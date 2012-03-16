# Initialization of simpatico

import _simpatico

from hoomd_script import globals

simulation = None

is_initialized = False;

def using_hoomd():
    global is_initialized
    global simulation

    print "Initialize Simpatico"
    simulation = _simpatico.HOOMDSimulation(globals.system_definition)
    is_initialized = True

