#ifndef SIMP_NOPAIR
#ifndef HOOMD_DPD_PAIR_H
#define HOOMD_DPD_PAIR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdPair.h"

#include <hoomd/PotentialPairGPU.h>
#include <hoomd/EvaluatorPairDPDThermo.h>
#include <hoomd/AllDriverPotentialPairGPU.cuh>

namespace McMd
{

   extern char classNameHoomdDpd[];

   /**
   * A potential encapsulating the HOOMD DPD evaluator
   *
   * \ingroup Pair_Module
   */
   class HoomdDpdPair : public HoomdPair< EvaluatorPairDPDThermo,
      gpu_compute_dpdthermo_forces, classNameHoomdDpd >
   {
   
   public:

      /**
      * Default constructor.
      */
      HoomdDpdPair();

      /**
      * Copy constructor
      */
      HoomdDpdPair(const HoomdDpdPair& other);

      /**
      * read parameters from file
      *
      * \param in input stream
      */
      void readParameters(std::istream &in);

      /**
      * Set DPD interaction energy for a specific pair of Atom types.
      *
      * \param i        type of Atom 1
      * \param j        type of Atom 2
      * \param epsilon  LJ energy parameter
      */
      void setEpsilon(int i, int j, double epsilon);

      /**
      * Get DPD interaction energy for a specific pair of Atom types.
      *
      * \param i   type of Atom 1
      * \param j   type of Atom 2
      * \return    epsilon_[i][j]
      */
      double epsilon(int i, int j) const;

      /**
      * Get DPD interaction cutoff distance squared for a specific pair of Atom types.
      *
      * \param i   type of Atom 1
      * \param j   type of Atom 2
      * \return    cutoffSq_[i][j]
      */
      double cutoffSq(int i, int j) const;

      /**
      * Modify a parameter, identified by a string.
      *
      * \param name   parameter name
      * \param i      atom type index 1
      * \param j      atom type index 2
      * \param value  new value of parameter
      */
      void set(std::string name, int i, int j, double value);

      /**
      * Get a parameter value, identified by a string.
      *
      * \param name   parameter name
      * \param i      atom type index 1
      * \param j      atom type index 2
      */
      double get(std::string name, int i, int j) const;

      /**
      * Get the Hoomd shift mode for this potential.
      */
      PotentialPairGPU<EvaluatorPairDPDThermo,
      gpu_compute_dpdthermo_forces>::energyShiftMode hoomdShiftMode() const
      {
         return PotentialPairGPU<EvaluatorPairDPDThermo,
            gpu_compute_dpdthermo_forces>::shift;
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
