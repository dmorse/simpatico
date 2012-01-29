// Include boost.python to do the exporting
#include <boost/python.hpp>
using namespace boost::python;

#include <mcMd/hoomdSimulation/HOOMDSimulation.h>
#include <mcMd/hoomdSimulation/HOOMDSystem.h>
#include <mcMd/potentials/pair/HOOMDAllPotentials.h>

using namespace McMd;

// specify the python module. Note that the name must expliclty match the PROJECT() name provided in CMakeLists
// (with an underscore in front)
BOOST_PYTHON_MODULE(_simpatico)
    {
#ifdef MCMD_HOOMD
    export_HOOMDSimulation();
    export_AllHOOMDPotentials();
#endif
    }

