#!/usr/bin/env python3

#    Copyright (C) 2017 Imperial College London and others.
#
#    Please see the AUTHORS file in the main source directory for a full list
#    of copyright holders.
#
#    Prof. C Pain
#    Applied Modelling and Computation Group
#    Department of Earth Science and Engineering
#    Imperial College London
#
#    amcgsoftware@imperial.ac.uk
#
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation,
#    version 2.1 of the License.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
#    USA

#Authors: C. Heaney, C.C. Pain, P. Salinas, and others.....

import sys
import libspud
import datetime
import csv

class Data_Assimilation():
    def __init__(self):
        self.fwd_input_file=""

class Halite_input():
    def __init__(self):
        self.output_filename = ""                       # fwd
        self.executable = ""                            # fwd
        self.input_file = ""                            # fwd 
        self.dimension = 0                              # fwd 
        self.avoid_perturbing_near_wells  = False       # Halite
        self.Halite_operation = ""                      # Halite
        self.functional = Functional()                  # func
        self.Field_To_study_list=[]
        ###   GA optimisation    ###
        self.ga_max_generations = 0
        self.ga_locations_to_study = 0
        self.ga_initial_guess = ""
        self.ga_gradient_convergence = -1
        self.ga_absolute_convergence = -1
        self.ga_population_generation = 0
        self.ga_breeding_prob = 0.0
        self.ga_mutation_prob = 0.0
        self.ga_fitness_functional = ""
        self.ga_Minimise = True
        self.ga_evol_algorithm = 1
        self.ga_selection_method = 1
        self.ga_mu = 0
        self.ga_lambda = 0
        self.ga_CMA = False
        self.ga_centroid = 0
        self.ga_sigma = 0
        self.ga_variables = []
        self.ga_hall_of_fame = 10
        self.number_processors = 1
        self.MPI_runs = 1
        self.precision = 0.01
        self.mutation_method = ""
        self.generation_method = ""
        self.mind_the_gap = -1
        self.cubit_path = ""
        self.cubit_input_file = ""
        self.cubit_pro_path = ""
        self.cubit_pro_input_file = ""
        self.producer_ids = ""
        self.injector_ids = ""
        self.optimise_input = False


class ga_variable():
    def __init__(self):
        self.name=""
        self.min_limit = 0
        self.max_limit = 0
        self.variable_pattern = ""
        self.normaliser = 1.0

def get_Halite_options():
    "Initialises the structure for the options."
    Halite_options = Halite_input()
    nirom_options = Nirom_input()
    fwd_options = Forward_model_options()

    #Read the input data from the user
    Halite_input_file = str(sys.argv[1])
    #Load the .Halite file
    libspud.load_options(Halite_input_file)
    ##Create Halite Variable
#    Halite_options = Halite_input()
    ##Start reading options
    print(("reading Halite input file", Halite_input_file))
    if  libspud.have_option('/avoid_perturbing_near_wells'):
        Halite_options.avoid_perturbing_near_wells = True

    
    Halite_options.Halite_operation = 'ga_optimisation'
    path = '/ga_optimisation'

    # previously these three options were at the top level
    Halite_options.output_filename = libspud.get_option(path + '/Output_filename')
    Halite_options.executable = libspud.get_option(path + '/Executable')
    Halite_options.input_file = libspud.get_option(path + '/Input_file')
    
    Halite_options.ga_max_generations = libspud.get_option(path+'/convergence_settings/Maximum_generations')

    Halite_options.ga_locations_to_study = libspud.get_option(path+'/Locations_to_study/value')
    if libspud.have_option(path+'/Locations_to_study/Initial_guess'):
        Halite_options.ga_initial_guess = libspud.get_option(path+'/Locations_to_study/Initial_guess')
    if libspud.have_option(path + '/Locations_to_study/Mind_the_gap'):
        Halite_options.mind_the_gap = libspud.get_option(path + '/Locations_to_study/Mind_the_gap')
    if libspud.have_option(path + '/cubit_integration'):
        Halite_options.cubit_path = libspud.get_option(path + '/cubit_integration/cubit_path')
        Halite_options.cubit_input_file = libspud.get_option(path + '/cubit_integration/cubit_input_file')
    if libspud.have_option(path + '/cubit_integration/cubit_pro'):
        Halite_options.cubit_pro_path = libspud.get_option(path + '/cubit_integration/cubit_pro/cubit_pro_path')
        Halite_options.cubit_pro_input_file = libspud.get_option(path + '/cubit_integration/cubit_pro/cubit_pro_input_file')            
    Halite_options.optimise_input = libspud.have_option(path + '/optimise_input')
    if libspud.have_option(path+'/convergence_settings/Gradient_convergence'):
        Halite_options.ga_gradient_convergence = libspud.get_option(path+'/convergence_settings/Gradient_convergence')
    if libspud.have_option(path+'/convergence_settings/Absolute_convergence'):
        Halite_options.ga_absolute_convergence = libspud.get_option(path+'/convergence_settings/Absolute_convergence')
    if libspud.have_option(path+'/Number_processors'):
        Halite_options.number_processors = libspud.get_option(path+'/Number_processors')
        if libspud.have_option(path + '/Number_processors/MPI_runs'):
            Halite_options.MPI_runs = libspud.get_option(path + '/Number_processors/MPI_runs')
            Halite_options.number_processors /= Halite_options.MPI_runs
    Halite_options.ga_population_generation = max(libspud.get_option(path+'/Population_generation'),2) #At least 2
    Halite_options.ga_breeding_prob = libspud.get_option(path+'/Breeding_probability/value')
    if libspud.have_option(path+'/Breeding_probability/Generation_method'):
        Halite_options.generation_method = libspud.get_option(path+'/Breeding_probability/Generation_method')
    Halite_options.ga_mutation_prob = libspud.get_option(path+'/Mutation_probability/value')
    if libspud.have_option(path+'/Mutation_probability/Mutation_method'):
        Halite_options.mutation_method = libspud.get_option(path+'/Mutation_probability/Mutation_method')
    Halite_options.precision = libspud.get_option(path + '/Precision')
    if libspud.have_option(path+'/Fitness/python_function'):
        Halite_options.ga_fitness_functional = libspud.get_option(path+'/Fitness/python_function')
    Halite_options.ga_Minimise = libspud.have_option(path+'/Fitness/Optimisation/Minimise')
    Halite_options.producer_ids = libspud.get_option(path+'/Fitness/producer_ids')
    Halite_options.injector_ids = libspud.get_option(path+'/Fitness/injector_ids')
    if libspud.have_option(path+'/GA_methods/Evolutionary_algorithm/eaSimple'):
        Halite_options.ga_evol_algorithm = 1
    elif libspud.have_option(path+'/GA_methods/Evolutionary_algorithm/eaMuCommaLambda'):
        Halite_options.ga_evol_algorithm = 2
        Halite_options.ga_mu = libspud.get_option(path+'/GA_methods/Evolutionary_algorithm/eaMuCommaLambda/Mu')
        Halite_options.ga_lambda = libspud.get_option(path+'/GA_methods/Evolutionary_algorithm/eaMuCommaLambda/Lambda')
    elif libspud.have_option(path+'/GA_methods/Evolutionary_algorithm/eaMuPlusLambda'):
        Halite_options.ga_evol_algorithm = 3
        Halite_options.ga_mu = libspud.get_option(path+'/GA_methods/Evolutionary_algorithm/eaMuPlusLambda/Mu')
        Halite_options.ga_lambda = libspud.get_option(path+'/GA_methods/Evolutionary_algorithm/eaMuPlusLambda/Lambda')
    if libspud.have_option(path+'/GA_methods/Selection_method/selBest'):
        Halite_options.ga_selection_method = 1
    elif libspud.have_option(path+'/GA_methods/Selection_method/selNSGA2'):
        Halite_options.ga_selection_method = 2
    elif libspud.have_option(path+'/GA_methods/Selection_method/selSPEA2'):
        Halite_options.ga_selection_method = 3

    Halite_options.ga_CMA = False
    if libspud.have_option(path+'/GA_methods/Use_CMA'):
        Halite_options.ga_CMA = True
        Halite_options.ga_centroid = libspud.get_option(path+'/GA_methods/Use_CMA/centroid')
        Halite_options.ga_sigma = libspud.get_option(path+'/GA_methods/Use_CMA/sigma')
    if libspud.have_option(path+'Hall_of_fame'):
        Halite_options.ga_hall_of_fame = libspud.get_option(path+'Hall_of_fame')
    nfields = libspud.option_count(path+'/Variables')
    for i in range(nfields):
        temp_field = ga_variable()
        fieldpath = path+'/Variables[' + str(i) + ']'
        temp_field.name = libspud.get_option(fieldpath+'/name')
        temp_field.min_limit = libspud.get_option(fieldpath+'/Min_limit')
        temp_field.max_limit = libspud.get_option(fieldpath+'/Max_limit')
        temp_field.variable_pattern = libspud.get_option(fieldpath+'/Variable_pattern')
        temp_field.normaliser = 1.0
        if libspud.have_option(fieldpath+'/normaliser'):
            temp_field.normaliser = libspud.get_option(fieldpath+'/normaliser')
        ##Now append to the variables list
        Halite_options.ga_variables.append(temp_field)

    libspud.clear_options()

    return Halite_options, fwd_options, nirom_options
