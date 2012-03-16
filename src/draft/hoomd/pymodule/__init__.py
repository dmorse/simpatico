# -*- coding: iso-8859-1 -*-
# this file exists to mark this directory as a python module

# need to import all submodules defined in this directory
import _simpatico

import init

from hoomd_script import globals
from hoomd_script import pair

def update_potentials():
    for compute in globals.forces:
        if isinstance(compute, pair.pair):
            name = 'HOOMD' + compute.cpp_force.__class__.__name__
            if name.endswith('GPU'):
                name=name.replace('GPU','')

            print name;
            if not init.is_initialized:
               print "\nError: Attempting to set up potentials when "
               print "      Simpatico is not yet initialized.\n"
               raise RuntimeError('Error updating potentials');

            cpp_potential = getattr(_simpatico,name)(init.simulation);
            print cpp_potential; 
            init.simulation.setPairPotential(cpp_potential);

def run():
    if not init.is_initialized:
        print "\nError: Attempting to run simulation but "
        print "      Simpatico is not yet initialized.\n"
        raise RuntimeError('Error updating potentials');

    update_potentials()
    init.simulation.syncFromHOOMD()
