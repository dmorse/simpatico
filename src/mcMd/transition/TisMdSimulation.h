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
      * Run equilibrium path for Tis algorithm.
      *
      * This function runs an MD simulation from the current state, and
      * terminates either after nStep MD steps or when the reaction
      * coordinate reaches the product region, whichever comes first.
      * The reaction coordinate is computed with the interval defined
      * by the ReactionAnalyzer.
      *
      * Upon return, nCrossing is the number of forward crossings of
      * the boundary of the reactant domain, i.e., the number of times
      * the reaction coordinate passes througt the reactantCoordinate
      * value from below, i.e., the number of sampled values that are
      * greater than the reactantCoordinate for which the previous
      * value was less than the reactantCoordinate.
      */
      void runTisEquilibrium(int nStep,
                             ReactionAnalyzer analyzer,
                             int& nCrossing);

      /**
      * Integrate a path sample for transition interface sampling.
      *
      * This function runs an MD simulation from the current system
      * state and terminates when the reaction coordinate defined 
      * by the ReactionAnalyzer at a sample becomes either less than 
      * the lowerCoordinate parameter or greater than upperCoordinate.
      * At each step at which the reaction coordinate is computed, 
      *
      */
      void runTisPath(double lowerCoordinate,
                      double upperCoordinate,
                      double sampleProbability,
                      ReactionAnalyzer& analyzer,
                      MdSnapShotArray& samples,
                      double& maxCoordinate,
                      bool& crossesLower);

      // Building blocks for Path Sampling

      /**
      * Reverse all velocities (define in system and/or Snapshot?)
      */ 
      void reverseSystemVelocities();

      // Shoot move - make a member function of SnapShot?

   };

}
#endif

