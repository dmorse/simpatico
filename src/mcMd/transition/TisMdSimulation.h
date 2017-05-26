#ifndef MCMD_TIS_MD_SIMULATION_H
#define MCMD_TIS_MD_SIMULATION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace McMd 
{

   class TisMdSimulation : public MdSimulation
   {

      /**
      * Integrate a path sample. (maybe define in subclass)
      */
      void runPath(MdSnapShot const & initial, 
                   bool reverse, 
                   MdPathAnalyzer& analyzer, 
                   MdTrajectory& trajectory);

      // Building blocks for Path Sampling

      /**
      * Add random velocities to a snapshot.
      */ 
      void reverseSystemVelocities();

   };
}
#endif

