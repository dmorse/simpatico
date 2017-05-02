#ifdef  MCMD_PERTURB
#ifdef SIMP_EXTERNAL
#ifndef SIMP_NOPAIR
#ifndef MCMD_MC_PAIR_EXTERNAL_PERTURBATION_H
#define MCMD_MC_PAIR_EXTERNAL_PERTURBATION_H


/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/perturb/LinearPerturbation.h>      // base class
#include <mcMd/neighbor/CellList.h>               // member
#include <mcMd/simulation/Simulation.h>
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/potentials/pair/McPairPotentialImpl.h>
#include <mcMd/potentials/external/ExternalPotentialImpl.h>
#include <mcMd/chemistry/Atom.h>
#include <util/ensembles/EnergyEnsemble.h>
#include <util/global.h>              
#include <util/containers/DMatrix.h>              // member template
#include <util/containers/DArray.h>               // member template
#include <util/containers/Pair.h>                 // member template parameter


namespace McMd
{

   using namespace Util;

   class McSystem;

   /**
   * A Perturbation in the pair interaction epsilon(0,1) for any pair
   * potential supporting setEpsilon() and in the external parameter for
   * any external potential supporting setExternalParameter().
   *
   * An McPairExternalPerturbation describes a Hamiltonian with a variable
   * pair interaction between species 0 and 1 and with a variable external
   * parameter. The perturbation parameter consists of the pair interaction 
   * parameter epsilon(0, 1) and the external parameter.
   *
   * \ingroup McMd_Perturb_Module
   */
   template < class PairInteraction, class ExternalInteraction >
   class McPairExternalPerturbation : public LinearPerturbation<McSystem>
   {

   public:

      /**
      * Constructor. 
      * 
      * \param system parent McSystem
      * \param size   number of systems (communicator size).
      * \param rank   id of this system (communicator rank).
      */
      McPairExternalPerturbation(McSystem& system, int size, int rank);

      /**
      * Destructor
      */
      virtual ~McPairExternalPerturbation();

      /**
      * Read parameter epsilon(0, 1) and external parameter from file.
      *
      * \param in input stream (file or std in).
      */ 
      virtual void readParameters(std::istream& in);
 
      /**
      * Set pair interaction parameter epsilon(0,1) and external parameter 
      * for this System.
      */
      virtual void setParameter();
      
      /**
      * Return the perturbation parameter for this System.
      *
      * \param i index of the perturbation parameter. 
      * In this case, i can either be 0 or 1 since nParameters = 2.
      */
      virtual double parameter(int i) const;

      /**
      * Return (0-1 pair energy) / ( kT *epsilon(0,1) ) if i is 0
      * or (external potential energy)/ ( kT *externalParameter ) if i is 1
      */
      virtual double derivative(int i) const;
      
      PairInteraction& pairInteraction() const;
      
      ExternalInteraction& externalInteraction() const;
      
   private:

      /// Neighbor array for internal use
      mutable CellList::NeighborArray neighbors_;

      /// Pointer to pair interaction
      mutable PairInteraction* pairInteractionPtr_;

      /// Pointer to external interaction
      mutable ExternalInteraction* externalInteractionPtr_;
      
      /*
      Number of perturbation parameters associated with a System.
      nParameters = 2 for McPairExternalPerturbation.
      */
      int nParameters_;
     
   };

   /// Implementation

   /*
   * Constructor.
   */
   template < class PairInteraction, class ExternalInteraction >
   McPairExternalPerturbation<PairInteraction, ExternalInteraction>::McPairExternalPerturbation(McSystem& system, int size, int rank)
    : LinearPerturbation<McSystem>(system, size, rank),
      pairInteractionPtr_(0),
      externalInteractionPtr_(0)
   {  setClassName("McPairExternalPerturbation"); }

   /*
   * Destructor.
   */
   template < class PairInteraction, class ExternalInteraction >
   McPairExternalPerturbation<PairInteraction, ExternalInteraction>::~McPairExternalPerturbation<PairInteraction, ExternalInteraction>()
   {}

   /*
   * Read epsilon(0,1) and external parameter from file
   */
   template < class PairInteraction, class ExternalInteraction >
   void McPairExternalPerturbation<PairInteraction, ExternalInteraction>::readParameters(std::istream& in)
   { 
      Perturbation::readParameters(in);
      nParameters_ = Perturbation::getNParameters();
   }

   /* 
   */
   template < class PairInteraction, class ExternalInteraction >
   PairInteraction& McPairExternalPerturbation<PairInteraction, ExternalInteraction>::pairInteraction() const
   {
      if (pairInteractionPtr_ == 0) {
         McPairPotential* pairPtr = &(system().pairPotential());
         McPairPotentialImpl< PairInteraction >* implPtr = 0;
         implPtr = dynamic_cast< McPairPotentialImpl< PairInteraction >* >(pairPtr);
         if (implPtr == 0) {
            UTIL_THROW("Failed dynamic cast of McPairPotential");
         }
         pairInteractionPtr_ = &implPtr->interaction();
      }
      return *pairInteractionPtr_;
   }
   
   /* 
   */
   template < class PairInteraction, class ExternalInteraction >
   ExternalInteraction& McPairExternalPerturbation<PairInteraction, ExternalInteraction>::externalInteraction() const
   {
      if (externalInteractionPtr_ == 0) {
         ExternalPotential* externalPtr = &(system().externalPotential());
         ExternalPotentialImpl< ExternalInteraction >* implPtr = 0;
         implPtr = dynamic_cast< ExternalPotentialImpl< ExternalInteraction >* >(externalPtr);
         if (implPtr == 0) {
            UTIL_THROW("Failed dynamic cast of ExternalPotential");
         }
         externalInteractionPtr_ = &implPtr->interaction();
      }
      return *externalInteractionPtr_;
   }


   /*
   * Set the parameter epsilon(0,1) and external parameter for this McSystem.
   */
   template < class PairInteraction, class ExternalInteraction >
   void McPairExternalPerturbation<PairInteraction, ExternalInteraction>::setParameter()
   { 
      pairInteraction().setEpsilon(0, 1, parameter_[0]);
      externalInteraction().setExternalParameter(parameter_[1]);
   }

   /* 
   * i = 0: Get the epsilon(0,1) from the parent System.
   * i = 1: Get the externalParameter from the parent System.
   */
   template < class PairInteraction, class ExternalInteraction >
   double McPairExternalPerturbation<PairInteraction, ExternalInteraction>::parameter(int i) const
   {  
      if (i >= nParameters_) {
         UTIL_THROW("perturbation parameter index is out of bounds");
      }
      double param = 0.0;
      if (i == 0) {
        param = pairInteraction().epsilon(0, 1); 
      } else if (i == 1) {
        param = externalInteraction().externalParameter(); 
      } 
      return param;
   }

   /*
   * i = 0: Returns pair energy for unlike pairs / (kT*epsilon(0,1))
   * i = 1: Returns external potential energy / ( kT *external parameter )
   */
   template < class PairInteraction, class ExternalInteraction >
   double McPairExternalPerturbation<PairInteraction, ExternalInteraction>::derivative(int i) const
   {  
      // Preconditions
      if (i >= nParameters_) {
         UTIL_THROW("perturbation parameter index is out of bounds");
      }
      double deriv = 0.0;
      if ( i == 0 ) {
         if (fabs(parameter_[i] - parameter(i)) > 1.0E-8) {
            UTIL_THROW("Pair perturbation parameter is not set correctly");
         }
         double pairEnergy, rsq;
         Atom  *jAtomPtr, *kAtomPtr;
         int    nNeighbor, nInCell;
         int    ic, nc, j, k, jId, kId, jType, kType;
         // Loop over cells
         pairEnergy = 0.0;
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
                        pairEnergy += pairInteraction().
                              energy(rsq, jAtomPtr->typeId(), kAtomPtr->typeId());

                     }

                  }

                } // secondary atoms

            } // primary atoms

         } // cells

         // Multiply by the temperature factor.
         if (system().energyEnsemble().isIsothermal()) {
             pairEnergy *= system().energyEnsemble().beta();
         } else {
            UTIL_THROW("Non isothermal ensemble for McPairPerturbation.");
         }

         deriv = pairEnergy/parameter_[i];
      } else if ( i == 1 ) {
         if (fabs(parameter_[i] - parameter(i)) > 1.0E-8) {
            UTIL_THROW("External perturbation parameter is not set correctly");
         }
         double externalEnergy;
         System::MoleculeIterator molIter;
         Molecule::AtomIterator atomIter;

         externalEnergy = 0.0;

         for (int iSpec=0; iSpec < system().simulation().nSpecies(); ++iSpec){
             for (system().begin(iSpec, molIter); molIter.notEnd(); ++molIter){
                  for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
                      externalEnergy += system().externalPotential().energy(atomIter->position(), atomIter->typeId());
                  }
             }
         }

         // Multiply by the temperature factor.
         if (system().energyEnsemble().isIsothermal()) {
            externalEnergy *= system().energyEnsemble().beta();
         } else {
            UTIL_THROW("Non isothermal ensemble");
         }

         deriv = -1.0*(externalEnergy/parameter_[i]);

      } 
      return deriv;

   }
   
}

#endif  
#endif  // #ifndef SIMP_NOPAIR
#endif  // #ifdef SIMP_EXTERNAL
#endif  // ifdef MCMD_PERTURB
