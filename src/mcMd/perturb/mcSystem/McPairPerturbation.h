#ifdef  MCMD_PERTURB
#ifndef MCMD_NOPAIR
#ifndef MC_PAIR_PERTURBATION_H
#define MC_PAIR_PERTURBATION_H


/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/perturb/LinearPerturbation.h>      // base class
#include <mcMd/neighbor/CellList.h>               // member

#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/potentials/pair/McPairPotentialImpl.h>
#include <mcMd/chemistry/Atom.h>
#include <mcMd/ensembles/EnergyEnsemble.h>

namespace McMd
{

   using namespace Util;

   class McSystem;

   /**
   * A Perturbation in the pair interaction epsilon(0,1) for any pair
   * potential supporting setEpsilon().
   *
   * An McPairPerturbation describes a Hamiltonian with a variable
   * pair interaction between species 0 and 1. The perturbation
   * parameter is the pair interaction parameter epsilon(0, 1).
   *
   * \ingroup Perturb_Module
   */
   template < class Evaluator >
   class McPairPerturbation : public LinearPerturbation<McSystem>
   {

   public:

      /**
      * Constructor. 
      */
      McPairPerturbation(McSystem& system);

      /**
      * Destructor
      */
      virtual ~McPairPerturbation();

      /**
      * Read parameter epsilon(0, 1) from file.
      *
      * \param in input stream (file or std in).
      */ 
      virtual void readParam(std::istream& in);
 
      /**
      * Set pair interaction parameter epsilon(0,1) for this System.
      */
      virtual void setParameter();

      /**
      * Return the pair parameter epsilon(0,1) for this System.
      * 
      * \param i index of the perturbation parameter. 
      * In this case, i is always 0 since nParameters = 1.
      */
      virtual double parameter(int i) const;

      /**
      * Return (0-1 pair energy) / ( kT *epsilon(0,1) )
      *
      * \param i index of the perturbation parameter. 
      * In this case, i is always 0 since nParameters = 1.
      */
      virtual double derivative(int i) const;

      /**
      * Return the pair potential evaluator.
      */
      Evaluator& evaluator() const;
  
   private:

      /// Neighbor array for internal use
      mutable CellList::NeighborArray neighbors_;

      /// Pointer to evaluator
      mutable Evaluator* evaluatorPtr_;

      /*
      Number of perturbation parameters associated with a System.
      nParameters = 1 for McPairPerturbation.
      */
      int nParameters_;

   };

   /// Implementation

   /*
   * Constructor.
   */
   template < class Evaluator >
   McPairPerturbation<Evaluator>::McPairPerturbation(McSystem& system)
    : LinearPerturbation<McSystem>(system),
      evaluatorPtr_(0)
   {}

   /*
   * Destructor.
   */
   template < class Evaluator >
   McPairPerturbation<Evaluator>::~McPairPerturbation<Evaluator>()
   {}

   /*
   * Read epsilon(0,1) from file
   */
   template < class Evaluator >
   void McPairPerturbation<Evaluator>::readParam(std::istream& in)
   {  
      Perturbation::readParameters(in); 
      nParameters_ = Perturbation::getNParameters();
   }

   /*
   */
   template < class Evaluator >
   Evaluator& McPairPerturbation<Evaluator>::evaluator() const
   {
      if (evaluatorPtr_ == 0) {
         McPairPotential* pairPtr = &(system().pairPotential());
         McPairPotentialImpl< Evaluator >* implPtr = 0;
         implPtr = dynamic_cast< McPairPotentialImpl< Evaluator >* >(pairPtr);
         if (implPtr == 0) {
            UTIL_THROW("Failed dynamic cast of McPairPotential");
         }
         evaluatorPtr_ = &implPtr->evaluator();
      }
      return *evaluatorPtr_;
   }

   /*
   * Set the parameter epsilon(0,1) for this McSystem.
   */
   template < class Evaluator >
   void McPairPerturbation<Evaluator>::setParameter()
   {  evaluator().setEpsilon(0, 1, parameter_[0]); }

   /* 
   * Get the tempering variable from the parent System.
   */
   template < class Evaluator >
   double McPairPerturbation<Evaluator>::parameter(int i) const
   {  
      if (i > nParameters_) {
         UTIL_THROW("perturbation parameter index is out of bounds");
      }
      return evaluator().epsilon(0, 1);
   }

   /*
   * Return pair energy for unlike pairs / (kT*epsilon(0,1))
   */
   template < class Evaluator >
   double McPairPerturbation<Evaluator>::derivative(int i) const
   {  
      // Preconditions
      if (i > nParameters_) {
         UTIL_THROW("perturbation parameter index is out of bounds");
      }
      if (fabs(parameter_[i] - parameter(i)) > 1.0E-8) {
         UTIL_THROW("Perturbation parameter is not set correctly");
      }
      if (!system().energyEnsemble().isIsothermal()) {
         UTIL_THROW("Non isothermal ensemble for McPairPerturbation.");
      }

      double energy, rsq;
      Atom  *jAtomPtr, *kAtomPtr;
      int    nNeighbor, nInCell;
      int    ic, nc, j, k, jId, kId, jType, kType;

      // Loop over cells
      energy = 0.0;
      nc = system().pairPotential().cellList().totCells(); 
      for (ic = 0; ic < nc; ++ic) {

         // Get array of neighbors_
         system().pairPotential().cellList().getCellNeighbors(ic, neighbors_, nInCell);
         nNeighbor = neighbors_.size();
  
         // Loop over primary atoms in this cell
         for (j = 0; j < nInCell; ++j) {
            jAtomPtr = neighbors_[j];
            jId      = jAtomPtr->id();
            jType    = jAtomPtr->typeId();
          
            // Loop over secondary atoms in this and neighboring cells
            for (k = 0; k < nNeighbor; ++k) {
               kAtomPtr = neighbors_[k];
               kId      = kAtomPtr->id();
               kType    = kAtomPtr->typeId();
     
               // Count each pair only once 
               if (kId > jId && jType != kType) {

                  // Exclude masked pairs
                  if (!jAtomPtr->mask().isMasked(*kAtomPtr)) {

                     rsq = system().boundary().
                         distanceSq(jAtomPtr->position(), kAtomPtr->position());

                     energy += evaluator().
                            energy(rsq, jAtomPtr->typeId(), kAtomPtr->typeId());

                  }

               }

            } // secondary atoms

         } // primary atoms

      } // cells

      // Multiply by the temperature factor.
      if (system().energyEnsemble().isIsothermal()) {
         energy *= system().energyEnsemble().beta();
      } else {
         UTIL_THROW("Non isothermal ensemble for McPairPerturbation.");
      }

      return energy/parameter_[i]; 
   }

}

#endif  
#endif  // #ifndef MCMD_NOPAIR
#endif  // ifdef MCMD_PERTURB
