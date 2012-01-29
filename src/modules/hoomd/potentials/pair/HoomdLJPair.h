#ifndef MCMD_NOPAIR
#ifndef HOOMD_LJ_PAIR_H
#define HOOMD_LJ_PAIR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdPair.h"

#include <hoomd/PotentialPair.h>
#include <hoomd/EvaluatorPairLJ.h>
#include <hoomd/AllDriverPotentialPairGPU.cuh>

namespace McMd
{

   extern char classNameHoomdLJ[];

   /**
   * A potential encapsulating the HOOMD LJ evaluator
   *
   * \ingroup Pair_Module
   */
   class HoomdLJPair : public HoomdPair< EvaluatorPairLJ,
      gpu_compute_ljtemp_forces, classNameHoomdLJ >
   {
   
   public:

      /**
      * Default constructor.
      */
      HoomdLJPair();

      /**
      * Copy constructor
      */
      HoomdLJPair(const HoomdLJPair& other);

      /**
      * read parameters from file
      *
      * \param in input stream
      */
      void readParam(std::istream &in);

      /**
      * Set LJ interaction energy for a specific pair of Atom types.
      *
      * \param i        type of Atom 1
      * \param j        type of Atom 2
      * \param epsilon  LJ energy parameter
      */
      void setEpsilon(int i, int j, double epsilon);

      /**
      * Get LJ interaction energy for a specific pair of Atom types.
      *
      * \param i   type of Atom 1
      * \param j   type of Atom 2
      * \return    epsilon_[i][j]
      */
      double epsilon(int i, int j) const;

      /**
      * Get the Hoomd shift mode for this potential.
      */
      PotentialPairGPU<EvaluatorPairLJ,
      gpu_compute_ljtemp_forces>::energyShiftMode hoomdShiftMode() const
      {
         return PotentialPairGPU<EvaluatorPairLJ,
            gpu_compute_ljtemp_forces>::shift;
      }

   private:

      /// Epsilon parameter
      double epsilon_[MaxAtomType][MaxAtomType];
      
      /// Sigma parameter
      double sigma_[MaxAtomType][MaxAtomType]; 
   };
  
}

#endif
#endif
