#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#Authors: P. Salinas
import operator

from Halite_classes import *
from deap import tools, base, creator, cma
from importance_map_tools import *
import time
import libspud
from exodus2gmsh import convertExodusII2MSH
from PIL import Image

aquifer_layers = 2

fp = '/home/hf116/Desktop/corrected_masked_imperial.jpg'
w_dim = 520 #width of imperial space
h_dim = 360 #height of imperial space
bounding_thickness = 500 #(so Halite_max and Halite_min values is based on this ie. min = 50 and max = 50 + w_dim)

#Read Halite options, for the time being just so we can call the same subroutines
# read in the Halite options
Halite_options = get_Halite_options()
# get xml extension (mpml or flml)
xml_extension = get_xml_extension(Halite_options.input_file)
# load the options for the forward model from xml file
libspud.load_options(Halite_options.input_file)
# Production file
prod_file = libspud.get_option('/simulation_name') + "_outfluxes.csv"
# Final time to extract from the csv file
final_time = libspud.get_option('/timestepping/finish_time')
#Retrieve extra variables required to run the optimisation
MIN_Halite = Halite_options.ga_variables[0].min_limit
MAX_Halite = Halite_options.ga_variables[0].max_limit
y_MIN_Halite = Halite_options.ga_variables[1].min_limit
y_MAX_Halite = Halite_options.ga_variables[1].max_limit
spatial_precision = int((abs(MIN_Halite) + abs(MAX_Halite)) * Halite_options.precision)
##Global variable to be used for gradient convergence
previous_convergence = [0.0, 0.0]
##Geothermal if we have temperature field
geothermal = libspud.have_option('material_phase['+str(0)+']/scalar_field::Temperature')
##Already visited studied
explored_locations = []

## TODO: 1) FIND OPTIMAL COMBINATION OF ALGORITHMS, FOR EXAMPLE DIFFERENT CROSSOVER?
## TODO: 1.5) IMPLEMENT SOME SORT OF HILL CLIMBER/GRADIENT DESCENT METHOD FOR CROSSOVER?


#Call the creators in the root of the program to ensure it works with scoop and in parallel
if  Halite_options.ga_Minimise:
    creator.create("FitnessMin", base.Fitness, weights=(-1,))
    creator.create("Individual", list, fitness=creator.FitnessMin)
else:
    creator.create("FitnessMax", base.Fitness, weights=(1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMax)

def analyse_results(instance):
    cwd = os.getcwd()
    foldername = get_folder_name(instance)
    try:
        timeList, prod, inject, prod_temp, inject_temp, TR_hotwell_quotient, TR_coldwell_quotient, TR_hotwell_divisor, TR_coldwell_divisor, hotwell_watts, coldwell_watts = get_productions_and_injections(prod_file, geothermal, foldername, cwd)
    except:
        os.chdir(cwd)
        print("WARNING: Failed to load production files.")
        timeList = 0.
        prod =0.
        inject=0.
        prod_temp=0.
        inject_temp =0.
        TR_hotwell_quotient =0.
        TR_coldwell_quotient =0.
        TR_hotwell_divisor =0.
        TR_coldwell_divisor =0.
        hotwell_watts=0.
        coldwell_watts=0.


    walltime= 0
    try:
        walltime = get_walltime(foldername, cwd)
    except:
        os.chdir(cwd)
        print("WARNING: Failed to obtain walltime")

    return timeList, prod, inject, prod_temp, inject_temp, walltime, TR_hotwell_quotient, TR_coldwell_quotient, TR_hotwell_divisor, TR_coldwell_divisor, hotwell_watts, coldwell_watts


def get_folder_name(instance):
    #To try to compact as much as possible the length of the folder while being unique
    #we convert the number representation to base64 from base10
    def get_digit(d):
        if 0 <= d <= 9:
            # 0 - 9
            c = 48 + d
        elif 10 <= d <= 35:
            # A - Z
            c = 55 + d
        elif 36 <= d <= 61:
            # a - z
            c = 61 + d
        elif d == 62:
            # -
            c = 45
        elif d == 63:
            # +
            c = 43
        return chr(c)

    def encode(n):
        out = []
        while n:
            n, r = n // 64, n % 64
            out.append(get_digit(r))
        while len(out) < 6:
            out.append('0')
        return ''.join(out)

    #To concatenate all the numbers we need them to be positive (no minus sign between them)
    ensure_positve = 0
    for var in Halite_options.ga_variables:
        ensure_positve = max(abs(var.min_limit), ensure_positve)    
    instance_bak = convert_instance(instance)
    foldername =  ''
    for val in instance_bak:
        foldername +=  str(int(float(val)+float(ensure_positve)))
    #We use base 64 representation of the numbers to try to reduce the length of the folder
    foldername = Halite_options.input_file[:-5]+"_"+ encode(int(float(foldername)))

    return foldername

def modify_experiment(input_file, foldername, instance, cwd):

    def modify_input_file(INPUT_FILE):
        # Read in the file
        with open(INPUT_FILE, 'r') as file:
            filedata = file.read()
        # Ensure that the path is correct
        filedata = filedata.replace( "PATH_TO_FOLDER" , os.getcwd() )
        
        # Convert to physical values
        instance_bak = convert_instance(instance)
        ##Substitute a given pattern by the instance value
        for i in range(Halite_options.ga_locations_to_study):
            for j in range(len(Halite_options.ga_variables)):
                normaliser = Halite_options.ga_variables[j].normaliser
                if j<2: #only convert coordinates
                    filedata = filedata.replace(Halite_options.ga_variables[j].variable_pattern + str(i + 1), str( float(instance_bak[j + len(Halite_options.ga_variables) * i])/float(normaliser) ))
                else:
                    filedata = filedata.replace(Halite_options.ga_variables[j].variable_pattern + str(i + 1), str(instance[j + len(Halite_options.ga_variables) * i]))

        # Overwrite file
        with open(INPUT_FILE, 'w') as file:
            file.write(filedata)


    ##Get into the folder
    os.chdir(foldername)
    #Modify cubit input file if requested
    cubit_success = False
    attempts = 0 
    max_wait = 120 #Since one wait is 60 seconds this is 2 hours
    cubit_pro_on =len(Halite_options.cubit_pro_input_file)
    while not cubit_success:
      #We try for two hours max
      if attempts > max_wait: exit 
      attempts += 1
      if len(input_file)>1:
          modify_input_file(input_file)
          #Run cubit to create new exodusII file given the instance /dev/null
          string = Halite_options.cubit_path.rstrip() + " -nographics -nojournal -batch " + Halite_options.cubit_input_file.rstrip() + " %s > ./cubit_log.txt"
          print("Creating model using cubit ...")
          os.system(string)
          
          if cubit_pro_on:
            string = Halite_options.cubit_pro_path.rstrip() + " -nographics -nojournal -batch " + Halite_options.cubit_pro_input_file.rstrip() + " %s > ./cubit_pro_log.txt"
            print("Creating mesh using cubit pro...")
            os.system(string)
          
          #Now convert the output .e file into .msh
          meshfile = libspud.get_option('/geometry/mesh::CoordinateMesh/from_file/file_name')
          #Check if cubit has been successful
          cubit_success = os.path.isfile(meshfile+".e")
          if (not cubit_success):
              print("############################################################################")
              print("ERROR: cubit failed to create a mesh for: ", foldername)
              if (attempts == 1):
                print("Just in case check that the template .jou file has the correct path for the step files, the PATH_TO_FILE string for the cd command and the variable patterns are implemented.")
              if (cubit_pro_on):#If we have a license it may fail to link, therefore we wait, otherwise we consider the error is different and we move on
                print("ERROR: waiting 60 seconds until trying again in case the floating license is being used.")
                print("Waited so far: ", str(attemps), " minutes...")
                time.sleep(60)
              else:
                exit()
              print("############################################################################")
          else:
            convertExodusII2MSH(meshfile)
            #Decompose the mesh is required
            if Halite_options.MPI_runs > 1: os.system("fldecomp -n " + str(Halite_options.MPI_runs) + " " + meshfile)

    #Modify the mpml file if requested
    if Halite_options.optimise_input: modify_input_file(Halite_options.input_file)


    ##Return to original path
    os.chdir(cwd)
    return

def get_walltime(foldername, cwd):
    from fluidity_tools import stat_parser as stat
    walltime =1e50
    ##Get into the folder
    os.chdir(foldername)
    output_name = libspud.get_option('/simulation_name')
    walltime = stat('./' + output_name + '.stat')["ElapsedWallTime"]["value"][-1]
    ##Return to original path
    os.chdir(cwd)
    return walltime

def get_productions_and_injections(prod_file, geothermal, foldername, cwd):
    import csv
    nPhases = libspud.option_count('/material_phase')
    ##Get into the folder
    os.chdir(foldername)
    String_id = "-S"
    String_prod = "- Volume rate"

    #Replace spaces by commas to make it consistent
    Halite_options.producer_ids.replace(' ', ',')
    Halite_options.injector_ids.replace(' ', ',')
    #Create a list with the integers
    producer_list = [int(e) if e.isdigit() else e for e in Halite_options.producer_ids.split(',')]
    injector_list = [int(e) if e.isdigit() else e for e in Halite_options.injector_ids.split(',')]

    with open(prod_file, 'r') as csvfile:
        datareader = csv.reader(csvfile, delimiter=',', quotechar='|')
        counter = -1 #To account for the header
        for row in datareader:
            try:
                counter+=1
            except:
                continue

    timeList = np.zeros(counter)
    prod = np.zeros((nPhases, counter))
    inject = np.zeros((nPhases, counter))
    prod_temp = np.zeros((nPhases, counter))
    inject_temp = np.zeros((nPhases, counter))
    counter = 0
    TR_hotwell_quotient = []
    TR_coldwell_quotient = []
    TR_hotwell_divisor = []
    TR_coldwell_divisor = []
    hotwell_watts = []
    coldwell_watts = []

    with open(prod_file, 'r') as csvfile:
        datareader = csv.reader(csvfile, delimiter=',', quotechar='|')
        prod_cols =[]
        inject_cols = []

        header = True
        previousT = 0.0
        for row in datareader:
            try:
                if header:
                    #Find the positions of interest of the producers
                    for ids in producer_list:
                        for i in range(len(row)):
                            if (String_id+str(ids)+String_prod) in row[i]:
                                prod_cols.append(i)
                    #Find the positions of interest of the injectors
                    for ids in injector_list:
                        for i in range(len(row)):
                            if (String_id+str(ids)+String_prod) in row[i]:
                                inject_cols.append(i)
                    header = False
                    continue

                phase_prod = 0
                phase_inject = 0
                timeList[counter] = float(row[0])
                # Calculate deltaT from absolute times
                deltaT = float(row[0]) - previousT
                previousT = float(row[0])

                current_time= float(row[0])
                t_mod = int(current_time/15768001)+(current_time%15768001<0)
                initial_producer_flag = bool(t_mod%2==0) #if flag is true then initial producer is producing
                #print('initial producer flag is', initial_producer_flag)

                #Initial temperature for hot well is 293.15
                #Initial temperature for cold well is 283.15
                #Let Tnatural = 287.15
                #T_hot - Tnatural = 6
                #T_cold - Tnatural = -4
                #heat capacity of water is 4185.5
            

                for i in range(len(row)):
                    # Update production
                    for j in range(len(prod_cols)):
                        if i == prod_cols[j]:
                            if geothermal:
                                #For geothermal I need to do temperature * production * rho * Cp
                                #CURRENTLY NOT INCLUDED Cp nor RHO
                                #First temperature, so we can calculate the production in this time level
                                prod_temp[phase_prod, counter] = abs(float(row[i]) * deltaT * float(row[i+nPhases*2]))
                            prod[phase_prod, counter] = abs(float(row[i]) * deltaT)

                            #prod is always surface 12
                            #when hotwell is producing
                            if (phase_prod == 1) and (initial_producer_flag == True) and (float(row[1]) >= 0.5):
                                #print('hot water produced with temperature difference', abs(287.15-float(row[i+nPhases*2]))) # and is', abs(float(row[i]) * deltaT * (293.15- float(row[i+nPhases*2]))))
                                #TR_hot_divisor = abs(float(row[i]) * deltaT *6)
                                TR_hotwell_quotient.append(abs(float(row[i]) * deltaT * (float(row[i+nPhases*2])-287.15)))
                                #print('divisor is', float(row[i]) * deltaT * 6.)
                                #print('here')

                                volume_rate = float(row[i])
                                mass_rate = volume_rate * 1000 
                                #calculate watts based on what is being produced
                                hotwell_watts.append(mass_rate * 4185.5 * (float(row[i+nPhases*2])-287.15))

                                #hotwell_watts.append( mass* 4185.5 * (float(row[i+nPhases*2])-287.15) )

                            #when hotwell is experiencing injection
                            if (phase_prod == 1) and (initial_producer_flag == False) and (float(row[1]) >= 0.5):
                                #print('temperature of water injected into hotwell is', float(row[i+nPhases*2]))
                                TR_hotwell_divisor.append(abs(float(row[i]) * deltaT * (float(row[i+nPhases*2])-287.15)))

                            phase_prod += 1
                            if phase_prod == nPhases: phase_prod = 0

                    #inject is always surface 11
                    # Update Injection
                    for j in range(len(inject_cols)):
                        if i == inject_cols[j]:
                            if geothermal:
                                #For geothermal I need to do temperature * production * rho * Cp
                                #CURRENTLY NOT INCLUDED Cp nor RHO
                                #First temperature, so we can calculate the production in this time level
                                inject_temp[phase_inject, counter] = abs(float(row[i]) * deltaT * float(row[i+nPhases*2]))

                            inject[phase_inject, counter] = abs(float(row[i]) * deltaT)

                            #when coldwell is producing
                            if (phase_inject == 1) and (initial_producer_flag == False) and (float(row[1]) >= 0.5):
                                #print('cold water produced with temperature difference', abs(287.15-float(row[i+nPhases*2]))) # and is', abs(float(row[i]) * deltaT * (283.15- float(row[i+nPhases*2]))))
                                #TR_coldwell.append(abs(float(row[i]) * deltaT * (283.15- float(row[i+nPhases*2]))))
                                TR_coldwell_quotient.append(abs(float(row[i]) * deltaT * (float(row[i+nPhases*2])-287.15)))

                                #print('temperature of water produced by coldwell is', float(row[i+nPhases*2]))

                                volume_rate = float(row[i])
                                mass_rate = volume_rate * 1000 
                                #calculate watts based on what is being produced
                                coldwell_watts.append(mass_rate * 4185.5 * abs(float(row[i+nPhases*2])-287.15))


                            #when coldwell is experiencing injection
                            if (phase_inject == 1) and (initial_producer_flag == True) and (float(row[1]) >= 0.5):
                                #print('temperature of water injected into coldwell is', float(row[i+nPhases*2]))
                                TR_coldwell_divisor.append(abs(float(row[i]) * deltaT * (float(row[i+nPhases*2])-287.15)))

                            phase_inject += 1
                            if phase_inject == nPhases: phase_inject = 0
                counter +=1
            except:
                continue
    ##Return to original path
    os.chdir(cwd)
    return timeList, prod, inject, prod_temp, inject_temp, TR_hotwell_quotient, TR_coldwell_quotient, TR_hotwell_divisor, TR_coldwell_divisor, hotwell_watts, coldwell_watts


#Function that evaluates the functional
#A penalty function can easily added here by providing a bad value if a requirement is not fulfilled
def fitness(instance):
    global  fitness_val
    timeList, prod, inject, prod_temp, inject_temp, walltime, TR_hotwell_quotient, TR_coldwell_quotient, TR_hotwell_divisor, TR_coldwell_divisor, hotwell_watts, coldwell_watts = analyse_results(instance)  # Only analyse
    if len(Halite_options.ga_fitness_functional) > 1:
        val = 0.
        try:
            #We need to define a global variable to connect the exce with the internal code
            exec("global  fitness_val;"+Halite_options.ga_fitness_functional)
            val = fitness_val
        except:
            print("#############################################################################")
            print("ERROR evaluating the user functional. Please check the functional introduced.")
            print(("Failure for:", get_folder_name(instance),". Value of 0 assigned."))
            print("#############################################################################")
    else:
        if geothermal:
            total_hotwell_TR_quotient = float(np.sum(TR_hotwell_quotient))
            total_coldwell_TR_quotient = float(np.sum(TR_coldwell_quotient))

            total_hotwell_TR_divisor = float(np.sum(TR_hotwell_divisor))

            # if total_hotwell_TR_divisor < 0.1:
            #     total_hotwell_TR_divisor = 0.1 
            #     print('time too short')

            total_coldwell_TR_divisor = float(np.sum(TR_coldwell_divisor))

            # if total_coldwell_TR_divisor < 0.1:
            #     total_coldwell_TR_divisor = 0.1 
            #     print('time too short')


            #average_hotwell_watts = float(float(np.sum(hotwell_watts))/len(hotwell_watts))

            #average_coldwell_watts = float(float(np.sum(coldwell_watts))/len(coldwell_watts))

            #print('len(coldwell_watts)', len(coldwell_watts))
            #print('coldwell_watts is',coldwell_watts)

            #print('total_coldwell_TR_quotient is',total_coldwell_TR_quotient)
            #print('total_coldwell_TR_divisor is',total_coldwell_TR_divisor)

            if total_hotwell_TR_divisor == 0.:
                
                cold_TR = total_coldwell_TR_quotient / total_coldwell_TR_divisor

                average_coldwell_watts = float(float(np.sum(coldwell_watts))/len(coldwell_watts))

                val = cold_TR * average_coldwell_watts #+ (hot_TR * total_coldwell_watts)

            elif total_coldwell_TR_divisor == 0.:
                hot_TR = total_hotwell_TR_quotient / total_hotwell_TR_divisor

                average_hotwell_watts = float(float(np.sum(hotwell_watts))/len(hotwell_watts))

                val = hot_TR * average_hotwell_watts
            
            else:
                cold_TR = total_coldwell_TR_quotient / total_coldwell_TR_divisor
                hot_TR = total_hotwell_TR_quotient / total_hotwell_TR_divisor

                average_hotwell_watts = float(float(np.sum(hotwell_watts))/len(hotwell_watts))

                average_coldwell_watts = float(float(np.sum(coldwell_watts))/len(coldwell_watts))

                val = ((cold_TR * average_coldwell_watts) + (hot_TR * average_hotwell_watts))/2
                

        
            
        else:
            val = np.sum(prod)
    return val,

def space_search_random(MINval, MAXval):

    val = (random.randint(MINval,MAXval)/spatial_precision ) * spatial_precision
    return val


from numpy import ones,vstack
from numpy.linalg import lstsq

def children(child1,child2):
    Nvar = len(Halite_options.ga_variables)
    Nwells = Halite_options.ga_locations_to_study
    new_child = []

    for i in range(Nwells):
        child1_well = child1[i * Nvar:i * Nvar + Nvar]
        child2_well = child2[i * Nvar:i * Nvar + Nvar]
  
        x_coords = (child1_well[0], child2_well[0])
        y_coords = (child1_well[1], child2_well[1])

        A = vstack([x_coords,ones(len(x_coords))]).T
        m, c = lstsq(A, y_coords)[0]

        random.seed(None)
        rand_x =random.uniform(min(child1_well[0],child2_well[0]),max(child1_well[0],child2_well[0]))
        new_child.append(round(rand_x))
        rand_y = m*rand_x + c 
        new_child.append(round(rand_y))

        print('new_child created', new_child)

        #rand_g_top = random.uniform(min(child1_well[2],child2_well[2]),max(child1_well[2],child2_well[2]))
        #new_child.append(round(rand_g_top))
        #rand_g_bottom = random.uniform(min(child1_well[3],child2_well[3]),max(child1_well[3],child2_well[3]))
        #new_child.append(round(rand_g_bottom))
    
        #conversion to binary format
        bin_well1 = format(int(child1_well[2]),"b")
        bin_well2 = format(int(child2_well[2]),"b")

        bin_well1 = bin_well1.zfill(64)
        bin_well2 = bin_well2.zfill(64)

        print('binwell1 is',bin_well1)
        print('binwell2 is',bin_well2)

        wells_open = False
        it_loop = 0

        while wells_open == False:

            bin_child = ''

            random.seed(None)
            selection_array = np.random.randint(0,2,64) #if 0 take from bin_well1 and if 1 take from bin_well2
            print('random num is', selection_array)

            for i in range(64):
                if selection_array[i] == 0:
                    bin_child = bin_child + str(list(bin_well1)[i])
                else:
                    bin_child = bin_child + str(list(bin_well2)[i])

            if it_loop==10000:
                for k in range(aquifer_layers):
                    list(bin_child)[-1*k] = 1
                        
            for j in range(aquifer_layers):
                    it_loop +=1
                    
                    print(it_loop)
                    if int(list(bin_child)[-1*j]) == 1:
                        bin_child = int(bin_child,2) #conversion into int
                        print('satisfied bin_child',bin_child)
                        new_child.append(bin_child)
                        wells_open = True
                      
    print('new child is in here is ',new_child)

    return new_child


def in_bounds(new_child,min_bound,max_bound):
    Nvar = len(Halite_options.ga_variables)
    Nwells = Halite_options.ga_locations_to_study

    num = 0

    print('new child is',new_child)

    for i in range(Nwells):
        one_well = new_child[i * Nvar:i * Nvar + Nvar]
        print('one well is', one_well)

        for j in range(Nvar):
            if j<2:
                if (min_bound<one_well[i]<max_bound):
                    num+=1
            else:
                if (0<new_child[i]<18446744073709551616):
                    num+=1

    if num == len(new_child):
        print('bounds satisfied')
        return True
    else:
        return False


def eaAlgorithm_by_steps(pop, toolbox, CXPB, MUTPB, NGEN, halloffame):


    print('pop is', pop)


    #We do not want to evaluate the first initial guesses as they will not be within the feasible domain

    # # Evaluate the entire population
    # counter = run_in_batches(pop)

    # # Evaluate
    # fitnesses = list(map(toolbox.evaluate, pop))
    # for ind, fit in zip(pop, fitnesses):
    #     ind.fitness.values = fit


    print('Halite_options.ga_population_generation is ', Halite_options.ga_population_generation)

    history_success_offspring = []

    for g in range(NGEN):

        if Halite_options.ga_CMA:
            #If CMA then generate new population like this
            offspring = toolbox.generate() 

            print('CMA is on')

        else: #Otherwise mate and mutation


            print('Else is on')

            # Select the next generation individuals
            if Halite_options.ga_evol_algorithm == 1:
                print('Equals to one with len', len(pop))
                offspring = toolbox.select(pop, len(pop))


            else:
                offspring = toolbox.select(pop, Halite_options.ga_lambda)
                print('Not one so lambda')
            # Clone the selected individuals
            offspring = list(map(toolbox.clone, offspring))
            print('offspring is', offspring)
            #Now we are applying a VarAnd method since we do, mate AND mutation
            # Apply crossover and mutation on the offspring

            success_offspring = toolbox.select(pop, 0)

            for child1, child2 in zip(offspring[::2], offspring[1::2]):

                feasibility_flag = 0

                r_count = 0

                while feasibility_flag < 2:
                    
                    #mate linearly
                    if r_count == 0: 
                        new_child = children(child1,child2)
                        new_child_copy = new_child.copy()
                    #expand search with radius dependent on r_count
                    else:
                        #new_child = children(child1,child2)
                        for i in range(2):
                            print('the random number is', random.uniform(-r_count, r_count))
                            new_child[i] += random.uniform(-r_count, r_count)
                            #new_child[i] += random.uniform(0, r_count)
                            new_child[i] = round(new_child[i])

                    #mutate (manual)
                    if random.random() <= CXPB:
                        index = random.randint(0, 1)
                        new_child[index] += space_search_random(MIN_Halite, MAX_Halite)
                        new_child[2] += random.randint(0, 18446744073709551615)

                    np_child = list(new_child)

                    print('feasibility is', feasible(np_child))

                    #dont need feasiblity check if sampling comes after mutation

                    if (feasible(np_child)) and (new_child not in history_success_offspring) and in_bounds(new_child,MIN_Halite,MAX_Halite):
                        print('this feasibility should be true:', feasible(np_child))
                        feasibility_flag+=1
                        print('converted success is', convert_instance(new_child))
                        deap_new_child = creator.Individual(new_child)
                        success_offspring.insert(0,deap_new_child)
                        history_success_offspring.append(new_child)
                        #restart r_count
                        r_count = 0

                    else:
                        r_count+=1
                        print('fail with r_count=', r_count)
                        #print('and failed child is', new_child)
                        #print('and child copy is', new_child_copy)
                        new_child = new_child_copy.copy()
                        time.sleep(5)

        print('success_offspring:', success_offspring)
        print('len(sucess_offspring) is', len(success_offspring))

        if Halite_options.ga_evol_algorithm == 2:
            pop[:] = toolbox.select(success_offspring, Halite_options.mu)
        elif Halite_options.ga_evol_algorithm == 3:
            pop[:] = toolbox.select(pop + success_offspring, Halite_options.mu)
        else:
            pop[:] = success_offspring

        if Halite_options.ga_CMA: toolbox.update(pop)  # Update strategy
        
        # #pop[:] = offspring
        # # The population is entirely replaced by the offspring
        # if Halite_options.ga_evol_algorithm == 2:
        #     pop[:] = toolbox.select(offspring, Halite_options.mu)
        # elif Halite_options.ga_evol_algorithm == 3:
        #     pop[:] = toolbox.select(pop + offspring, Halite_options.mu)
        # else:
        #     pop[:] = offspring

        # if Halite_options.ga_CMA: toolbox.update(pop)  # Update strategy


        print('NGEN is', NGEN)

        #Halite_options.ga_population_generation = len(pop)
        #Halite_options.number_processors = int(len(pop)/2)
        #Halite_options.ga_population_generation = int(ceil(float(Halite_options.ga_population_generation)/float(parProcess))  * parProcess)        

        #print('after feasible, len of pop is', len(pop))

        #Evaluate the individuals with an invalid fitness
        counter = run_in_batches(pop)

        ############This section seems to try to evaluate the new individuals to early##############
        #invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        #fitnesses = map(toolbox.evaluate, invalid_ind)
        #for ind, fit in zip(invalid_ind, fitnesses):
        #    ind.fitness.values = fit

        #Update hall of fame
        halloffame.update(pop)
        best_so_far = halloffame[0]
        best_distance_so_far = fitness(best_so_far)
        # Check current state of convergence
        #terminate(pop,best_distance_so_far, g)
        #Print best result so far:
        best_so_far_bak = convert_instance(best_so_far)
        print(("Best simulation at generation " + str(g+1) + ": ", best_so_far_bak, "; Produces: ", best_distance_so_far[0]))
        #Update as well best hall of fame
        create_output_file(halloffame)
        #create_best_folders(halloffame)

def One_point_crossover(ind1, ind2):
    size = min(len(ind1), len(ind2))
    if size == 1:
        point =1
    else:
        point = random.randint(1, size - 1) #if size == 1 this gives serious problems
    ind1[point:], ind2[point:] = ind2[point:], ind1[point:]
    return ind1, ind2


def run_in_batches(sub_pop):
    from subprocess import Popen
    from copy import deepcopy
    import time

    #For the time being, do not rerun simulations
    sub_population = deepcopy(sub_pop)
    pop_to_run = []
    processes = []
    active_run = []
    #Check that the locations haven't been checked already
    # for instance in sub_population:
    counter = 0
    for instance in sub_population:
        foldername = get_folder_name(instance)
        #We check that a simulation has been run if there is a production file already in place
        if os.path.isfile(foldername+"/"+prod_file):
            explored_locations.append(instance)

        if instance in explored_locations:
            try:
                pop_to_run.remove(instance)
            except:
                pass
        else:
            counter += 1
            pop_to_run.append(instance)
            explored_locations.append(instance)
        #Once we fill up the total number of CPUs, we continue
        #if len(pop_to_run) == Halite_options.number_processors: break

    
    #Create folders
    folders = []
    input_file = Halite_options.cubit_input_file
    # if libspud.have_option('/Halite_operation/ga_optimisation/cubit_integration'):
    #     input_file = Halite_options.cubit_input_file#Halite_options.input_file[:-5]
    cwd = os.getcwd()
    #Generate the sub-folders and create all the input files
    for instance in pop_to_run:
        cwd = os.getcwd()
        have_wells = True
        foldername = get_folder_name(instance)
        folders.append(foldername)
        existed = create_folder_to_run(Halite_options, cwd, have_wells, foldername)
        if existed:#When rerunning simulations with folders already in place
            sub_population.remove(instance)
            continue
        modify_experiment(input_file, foldername, instance, cwd)


    if len(pop_to_run)>=1 and len(folders)>0:
        # Prepare the commands to run the batch of simulations
        commands = []
        # Run the simulations in parallel
        for k in range(len(folders)):
            if Halite_options.MPI_runs > 1:
                string = "mpirun -n "+str(Halite_options.MPI_runs)+ " " + Halite_options.executable + " " + Halite_options.input_file
            else:
                string = Halite_options.executable + " " + Halite_options.input_file
            commands.append('cd ' + folders[k] + ' && ' + string)

        #Now time to run simulations
        for _ in range(len(commands)):
            # Check for finished jobs
            while not len(processes) < Halite_options.number_processors:
                time.sleep(60)
                for i, p in enumerate(processes):
                    if p.poll() != None:
                        print(('Run '+str(active_run[i])+' finished! status:' + str(p.poll())))
                        try:
                            processes[i].kill()
                        except:
                            pass
                        #Now remove the same info from all the lists to ensure that consistency is kept
                        processes.pop(i)
                        active_run.pop(i)
                        break
            # Run case 
            i = len(commands)-1
            processes += [subprocess.Popen(cmd, shell=True) for cmd in commands[i:i+1]]
            active_run.append(folders)
            commands.pop(i)
            folders.pop(i)

        # wait for completion for the remaining cases
        for p in processes: p.wait()
        #Cleanup the list
        del commands[:]
        
    #Just in case return to the root folder
    os.chdir(cwd)

    return counter

#Mutate method for integer
def mutate(instance, mutpb):

    if len(Halite_options.mutation_method)>1:
        #Execute python code from the user
        exec(Halite_options.mutation_method)
    else:
        #Wether to do the mutation or not is decided outside
        if random.random() <= 1.0:
            index = random.randint(0, len(instance) - 1)
            instance[index] += space_search_random(MIN_Halite, MAX_Halite)


    return instance,


# def checkBounds(MIN_Halite, MAX_Halite):
#     #Ensure that the results are bounded
#     # If many locations are studied and a clear distance is specified, then ensure that
#     # this is satisfied by the locations
#     Nvar = len(Halite_options.ga_variables)
#     Nwells = Halite_options.ga_locations_to_study
#     def decorator(func):
#         def wrapper(*args, **kargs):

#             offspring = func(*args, **kargs)
#             for child in offspring:
#                 #Ensure that wells are not too close
#                 # Perform this by "pushing" wells that are too close
#                 if Halite_options.mind_the_gap > 0:
#                     # First two variables of a location to study are the X and Y locations
#                     for i in range(Nwells - 1):
#                         Xorig = np.asanyarray(child[i * Nvar:i * Nvar + Nvar])
#                         Xorig = Xorig[0:2] #take only coordinates 
#                         for j in range(Nwells):
#                             if j == i: continue  # Ignore itself!
#                             Xother = np.asanyarray(child[j * Nvar:j * Nvar + Nvar])
#                             Xother = Xother[0:2] #take only coordinates
#                             dist = (np.dot(Xorig - Xother, Xorig - Xother)) ** 0.5
#                             if dist < Halite_options.mind_the_gap:
#                                 if abs(dist) < spatial_precision:  # Move the node away
#                                     Xother[0] = Xother[0] + (-1) ** random.randrange(2) * Halite_options.mind_the_gap
#                                     Xother[1] = Xother[1] + (-1) ** random.randrange(2) * Halite_options.mind_the_gap
#                                     #Xother[2] = Xother[2] + (-1) ** random.randrange(2) * Halite_options.mind_the_gap
#                                     #Xother[3] = Xother[3] + (-1) ** random.randrange(2) * Halite_options.mind_the_gap
#                                 else:
#                                     # Move node to new position
#                                     Xother = (Xother - Xorig) * Halite_options.mind_the_gap / dist + Xother
#                                 # Ensure that the node is in the space of search
#                                 print('Xother0 is', Xother[0])
#                                 print('Xother1 is', Xother[1])
#                                 print('j * Nwells is', j * Nwells)
#                                 #child[j * Nwells] = (np.int(Xother[0]) / spatial_precision) * spatial_precision
#                                 #child[j * Nwells + 1] = (np.int(Xother[1]) / spatial_precision) * spatial_precision
                
                
#                 for i in range(Nwells):
#                     #for x coord
#                     if child[i * Nvar] > MAX_Halite:
#                         child[i * Nvar] = MAX_Halite
#                     if child[i * Nvar] > MAX_Halite:
#                         child[i * Nvar] = MAX_Halite
                
#                 for i in range(Nwells):
#                     #for y coord
#                     if child[i * Nvar + 1] < MIN_Halite:
#                         child[i * Nvar + 1] = MIN_Halite
#                     if child[i * Nvar + 1] < MIN_Halite:
#                         child[i * Nvar + 1] = MIN_Halite
                    
#                 for i in range(Nwells):
#                     #for gamma
#                     if child[i * Nvar + 2] < 1:
#                         child[i * Nvar + 2] = 1
#                     elif child[i * Nvar + 2] > 18446744073709551615:
#                         child[i * Nvar + 2] = 18446744073709551615
#             return offspring
#         return wrapper
#     return decorator



def convert_instance(instance):
    #j =0
    instance_back = copy.deepcopy(instance)

    #Nvar = len(Halite_options.ga_variables)
    Nwells = Halite_options.ga_locations_to_study

    for j in range(Nwells): # convert only coordinates
        instance_back[j*Nwells] = linear_converter(instance_back[j*Nwells], MIN_Halite, MAX_Halite, Halite_options.ga_variables[0].min_limit, Halite_options.ga_variables[0].max_limit)
        instance_back[j*Nwells +1] = linear_converter(instance_back[j*Nwells+1], MIN_Halite, MAX_Halite, Halite_options.ga_variables[1].min_limit, Halite_options.ga_variables[1].max_limit)

    # for i in range(len(instance)):
    #     print('instance initially is', instance_back[i], 'with i:',i , 'max limit=', Halite_options.ga_variables[j].max_limit)
    #     #if i == 0: continue #Ignore first because it is the reference
    #     instance_back[i] = linear_converter(instance_back[i], MIN_Halite, MAX_Halite,
    #                      Halite_options.ga_variables[j].min_limit, Halite_options.ga_variables[j].max_limit)
    #     j += 1
    #     #Restart j so it iterates over values per location
    #     print('now instance back is', instance_back[i])
    #     if j%(len(instance_back)/ Halite_options.ga_locations_to_study)==0: j = 0
    
    return instance_back

def convert_instance_back(instance):
    #j = 0
    instance_back = copy.deepcopy(instance)

    Nwells = Halite_options.ga_locations_to_study

    for j in range(Nwells): # convert only coordinates
        instance_back[j*Nwells] = linear_converter(instance_back[j*Nwells], Halite_options.ga_variables[0].min_limit,Halite_options.ga_variables[0].max_limit, MIN_Halite, MAX_Halite)
        instance_back[j*Nwells +1] = linear_converter(instance_back[j*Nwells+1], Halite_options.ga_variables[1].min_limit,Halite_options.ga_variables[1].max_limit, MIN_Halite, MAX_Halite)

    # for i in range(len(instance)):
    #     #if i == 0: continue #Ignore first because it is the reference
    #     instance_back[i] = linear_converter(instance_back[i], Halite_options.ga_variables[j].min_limit,Halite_options.ga_variables[j].max_limit, MIN_Halite, MAX_Halite)
    #     j += 1
    #     #Restart j so it iterates over values per location
    #     if j%(len(instance_back)/ Halite_options.ga_locations_to_study)==0: j = 0
    
    return instance_back



def linear_converter(old_var, old_min, old_max, new_min, new_max):
    if new_max == old_max and new_min == old_min:
        val = old_var
    else:
        val = int((float(old_var - old_min) / float(old_max - old_min) * float(new_max - new_min) + new_min))
    return val

def get_converted_grid(image_fp, w_dim, h_dim, bounding_thickness):
    image = Image.open(image_fp)
    newsize = (520, 360)
    image = image.resize(newsize)
    print('size is', image.size)
    #choose a threshold value that will split image into built up areas and non-built up areas
    thresh = 40
    fn = lambda x : 0 if x > thresh else 255
    #convert image to greyscale then change it into binary image using a threshold value
    r = image.convert('L').point(fn)
    r = r.rotate(180).transpose(Image.FLIP_LEFT_RIGHT)
    r = r.resize((w_dim,h_dim))
    r_array = np.array(r)

    #space the will hold the converted coordinates
    bin_grid = np.zeros((w_dim,w_dim), dtype=bool)

    #get the spatial precision in the x-space to one in the y-space (step size)
    first_converted = linear_converter(int(MAX_Halite/2), MIN_Halite, MAX_Halite, y_MIN_Halite, y_MAX_Halite)
    second_converted = linear_converter(int(MAX_Halite/2)+spatial_precision, MIN_Halite, MAX_Halite, y_MIN_Halite, y_MAX_Halite)
    step_size = second_converted - first_converted
    print('step size is', step_size)
    print('spatial precision:', spatial_precision)

    i = 0
    j = 0

    for y in range(0,h_dim,step_size):
      for x in range(0,w_dim):
        if r_array[i*step_size][j] == 0:
          bin_grid[i*spatial_precision][j] = True
        j+=1
      i+=1
      j=0
    
    #add padding to the space - NEEDS TO BE COMPATIBLE WITH THE CUBIT MODEL ie. 
    # h_dim + 2* bounding_thickness = y (height) of model
    bin_grid = np.pad(bin_grid, bounding_thickness, mode='constant')

    return bin_grid

conv_grid = get_converted_grid(fp,w_dim,h_dim,bounding_thickness)

coords = np.argwhere(conv_grid==True)

coords_copy = coords.copy()

# y values
coords_copy[:,1] = coords[:,0] 

# x values
coords_copy[:,0] = coords[:,1] 

import math
from numpy.random import choice
coords_list = list(range(0, len(coords)))

#function to create child
#probability distribution based on distance (for single parent that has one child)
def create_child_with_probability(parent):
  Nvar = len(Halite_options.ga_variables)
  Nwells = Halite_options.ga_locations_to_study
  
  new_child = []

  for i in range(Nwells):
    well = parent[i * Nvar:i * Nvar + Nvar]
  
    distance_array = []
    e_array = []
  
    for n in range(len(coords)):
      x_dis = (well[0] - coords[n,1])**2
      y_dis = (well[1] - coords[n,0])**2
      distance = (x_dis+y_dis)**0.5
      distance_array.append(distance)
      e_array.append((10**(-distance)) * 1e250)

    sum_e = np.sum(e_array)
    prob_array = np.true_divide(e_array, sum_e)
    draw = choice(coords_list, 1, p=prob_array)
    
    for p in range(Nvar-1,-1,-1):
      print(coords[draw][0][p])
      new_child.append(coords[draw][0][p]) #,1],coords[draw,0]

  return new_child

#for depth-included wells, need to check order of coords!!!
# def old_feasible(instance):
#     Nvar = len(Halite_options.ga_variables)
#     Nwells = Halite_options.ga_locations_to_study
#     overall_N_count = 0
#     for i in range(Nwells):
#         N_count = 0
#         well = instance[i*Nvar: i*Nvar + Nvar]
#         if well[0] in coords_copy[:,0]:
#             index_values = np.argwhere(coords_copy[:,0] == well[0])
#             for n in range(spatial_precision,spatial_precision): #could use int(spatial_precision/2) instead
#               if well[1]+n in coords_copy[:,1][index_values]:
#                 N_count+=1
        
#         if N_count> 0:
#           overall_N_count+=1

#     if overall_N_count == Nwells:
#         return True

#     else:
#         return False

def feasible(instance):
    Nvar = len(Halite_options.ga_variables) #3
    Nwells = Halite_options.ga_locations_to_study #2
    overall_N_count = 0
    for i in range(Nwells):
        N_count = 0
        well = instance[i*Nvar: i*Nvar + Nvar]
        if well[0] in coords_copy[:,0]:
            #print(well[0])
            index_values = np.argwhere(coords_copy[:,0] == well[0])
            #print(coords_copy[:,1][index_values])
            for n in range(-spatial_precision,spatial_precision):
              if well[1]+n in coords_copy[:,1][index_values]: ##spational precision check? (np.int(Xother[0]) / spatial_precision) * spatial_precision
                N_count+=1
        
        if N_count> 0:
          overall_N_count+=1
            # for j in range(1,Nvar):
            #     if well[j] in coords_copy[:,j][index_values]:
            #         N_count+=1 
        
        #print(overall_N_count)

    if overall_N_count == Nwells:
        return True

    else:
        return False

# def feasible(instance):
#     return True

# Drop the first element,
# which will be replace by our initial guess.
def set_initial_guess(population, init_guess_string):
    #Convert string from diamond into list with integers
    initial_guess = [int(e) if e.isdigit() else e for e in init_guess_string.split(',')]
    if len(initial_guess) != len(Halite_options.ga_variables)*Halite_options.ga_locations_to_study:
        print("ERROR: The initial guess must be of the size of locations to study times number of variables")
        print("Initial guess ignored.")
        return
    population.pop()
    #Convert to internal space
    initial_guess = convert_instance_back(initial_guess)
    guess_ind = creator.Individual(initial_guess)
    #Make sure it is in the correct space

    population.insert(0, guess_ind)

####################CHANGE CONVERGENCE CRITERIA######################
#Finds the fittest
def get_best_result(population):

    if isinstance(population[0], list):
        fitness_values = list(map(fitness, population))
        if Halite_options.ga_Minimise:
            index = fitness_values.index(min(fitness_values))
        else:
            index = fitness_values.index(max(fitness_values))
        return population[index]
    else:
        if Halite_options.ga_Minimise:
            return min(population, key=operator.attrgetter('fitness'))
        else:
            return max(population, key=operator.attrgetter('fitness'))



def terminate(population, best_distance_so_far, generation):
    global previous_convergence
    reslt =fitness(get_best_result(population))
    #reslt = float('.'.join(str(ele) for ele in reslt1))
    if Halite_options.ga_Minimise:
        result = min(reslt[0],best_distance_so_far[0])
    else:
        result = max(reslt[0],best_distance_so_far[0])
    if Halite_options.ga_gradient_convergence > 0.:
        if abs(2.*result - sum(previous_convergence[:]))/max(abs(result),1e-16) < Halite_options.ga_gradient_convergence:
            if generation>1: raise StopIteration
        else:
            previous_convergence[1] = previous_convergence[0]
            previous_convergence[0] = result

    if Halite_options.ga_absolute_convergence > 0.:
        if Halite_options.ga_Minimise and result <= Halite_options.ga_absolute_convergence:
            raise StopIteration
        elif not Halite_options.ga_Minimise and result >= Halite_options.ga_absolute_convergence:
            raise StopIteration
    return False

#I presume it evaluates some sort of residual
def distance_from_best_result(population):
    result = get_best_result(population)
    return fitness(result)[0]

#####################################################################


#Outputs the best
def output(best_instance):
    best_instance_bak = convert_instance(best_instance)
    print(('Best result:', best_instance_bak))
    distance = fitness(best_instance)
    if distance[0] > 0.1:
        print(("The best result produces: ", abs(distance[0])))

#This functions creates an output file in csv format that includes the information of the best results
#Untested
def create_output_file(halloffame):
    import csv

    HallOfFame_file = "Best_results_"+Halite_options.output_filename+".csv"

    with open(HallOfFame_file, 'w') as csvfile:
        datarwriter = csv.writer(csvfile, delimiter=',', quotechar='|')
        #First write header
        writeList = ["Ranking"]
        writeList.append("Fitness evaluation")
        for k in range(Halite_options.ga_locations_to_study):
            for variable in Halite_options.ga_variables:
                writeList.append("Location " + str(k+1)+ ": " + variable.name)

        writeList.append("Folder name")
        datarwriter.writerow(writeList)
        #Now proceed to write the data
        i = 0
        for instance in halloffame:
            i += 1
            instance_bak = convert_instance(instance)
            writeList = [i] + [fitness(instance)[0]] + instance_bak + [get_folder_name(instance)]
            datarwriter.writerow(writeList)

    return

def create_best_folders(halloffame):


    #Create a folder with the best result
    try:
        foldername = get_folder_name(halloffame[0])
    except:
        foldername = ""
    #New name for the folder
    newfoldername = "Best_result"
    os.system("cp -r "+ foldername + " "+ newfoldername)

    vtufile = libspud.get_option('/simulation_name')+"_1.vtu"
    #Now we copy the initial vtu of all the results in the hall of fame, to be able to easily explore
    # the different well locations as well as the production .csv files
    newfoldername = "Best_well_configurations"
    os.system("mkdir " + newfoldername)
    for k in range(len(halloffame)):
        foldername = get_folder_name(halloffame[k])
        #Copy vtu file in order
        os.system("cp "+ foldername + "/" +vtufile + " " + newfoldername+ "/" +"best_result_ranked_" + str(k+1)+".vtu")
        #Copy production .csv fill
        os.system("cp "+ foldername + "/" +prod_file + " " + newfoldername+ "/" +"best_result_ranked_" + str(k+1)+".csv")

    return

#Initialise the genetic algorithm
def setup():

    toolbox = base.Toolbox()
    toolbox.register("attribute", space_search_random, MIN_Halite, MAX_Halite)
    toolbox.register("individual", tools.initRepeat, creator.Individual,
                     toolbox.attribute, n=len(Halite_options.ga_variables)*Halite_options.ga_locations_to_study)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    #toolbox.register("mate", tools.cxOnePoint)
    #toolbox.register("mutate", mutate, mutpb=Halite_options.ga_mutation_prob)
    #toolbox.decorate("mate", checkBounds(MIN_Halite, MAX_Halite))
    #toolbox.decorate("mutate", checkBounds(MIN_Halite, MAX_Halite))


    if Halite_options.ga_selection_method == 1:
        toolbox.register("select", tools.selBest)
    elif Halite_options.ga_selection_method == 2:
        toolbox.register("select", tools.selNSGA2)
    else:
        Halite_options.toolbox.register("select", tools.selSPEA2)

    toolbox.register("evaluate", fitness)

    #To parallelise using SCOOP
    #toolbox.register("map", futures.map)

    if Halite_options.ga_CMA:
        #Covariance Matrix Adaptation Evolution Strategy (CMA-ES) [Hansen2001]
        strategy = cma.Strategy(centroid=Halite_options.ga_centroid, sigma=Halite_options.ga_sigma)
        toolbox.register("generate", strategy.generate, creator.Individual)
        toolbox.register("update", strategy.update)


    return toolbox


def main():


    #The population of a generation has to be a multiple of the number of cpus used,  not doing this is un-optimal
    parProcess = Halite_options.number_processors
    Halite_options.ga_population_generation = int(ceil(float(Halite_options.ga_population_generation)/float(parProcess))  * parProcess)


    if not os.path.isfile(Halite_options.executable):
        print("#############################################################################")
        print("ERROR: IC-FERST executble not found in the given path.")
        print("#############################################################################")
        exit()
    
    toolbox = setup()
    #Create an initial population
    population = toolbox.population(n=Halite_options.ga_population_generation)
    #population = toolbox.population(n=100)
    #Specify that an initial seed value
    if Halite_options.ga_initial_guess:
        set_initial_guess(population, Halite_options.ga_initial_guess)
    stats = tools.Statistics()
    stats.register("best_instance_of_population", get_best_result)
    stats.register("distance", distance_from_best_result)
    stats.register("terminate", terminate)
    #halloffame = tools.HallOfFame(Halite_options.ga_population_generation*Halite_options.ga_max_generations)
    halloffame = tools.HallOfFame(min(Halite_options.ga_hall_of_fame,
                                    Halite_options.ga_population_generation* Halite_options.ga_max_generations-1))
    try:

        eaAlgorithm_by_steps(population, toolbox, Halite_options.ga_breeding_prob, Halite_options.ga_mutation_prob,
                                       Halite_options.ga_max_generations, halloffame)

    except StopIteration:
       print('Sad')
       pass
    finally:
        #Create list with best results
        create_output_file(halloffame)
        create_best_folders(halloffame)
        best_instance = halloffame[0]
        output(best_instance)
        return best_instance


# if __name__ == '__main__':
#     main()