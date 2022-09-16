# How-to:

To call the program you need python3 and call opal.py with an input file *.opal 
For example: python3 opal.py test.opal

To run the software it is imperature that the spud library contained in the folder spud is compiled and the bindings for python3 created. For this the following commands need to be run:
1) Open a terminal inside spud folder
2) Run: sudo ./configure
3) Run: sudo make 
4) Run: cd ./python
5) Run: sudo python3 setup.py build 
6) Run: sudo python3 setup.py install

This will create the necessary library for Opal to work properly.
If still libspud fails to be loaded, make sure that fluidity-dev and python3-spud libraries are installed properly. Another option is to perform sudo make install in the root of the spud folder.






