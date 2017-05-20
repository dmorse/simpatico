#ifdef  MCMD_PERTURB
#ifdef SIMP_EXTERNAL

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McExternalPerturbation.h"
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/chemistry/Atom.h>
#include <util/ensembles/EnergyEnsemble.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   template < class Interaction >
   McExternalPerturbation<Interaction>::McExternalPerturbation(McSystem& system,
                                                             int size, int rank)
    : LinearPerturbation<McSystem>(system, int size, int rank),
      interactionPtr_(0)
   { }

   /*
   * Destructor.
   */
   template < class Interaction >
   McExternalPerturbation<Interaction>::~McExternalPerturbation<Interaction>()
   { }

   /*
   * Read external parameter from file
   */
   template < class Interaction >
   void McExternalPerturbation<Interaction>::readParameters(std::istream& in)
   {  
      Perturbation::readParameters(in); 
      nParameters_ = Perturbation::getNParameters();
   }

   /*
   */
   template < class Interaction >
   Interaction& McExternalPerturbation<Interaction>::interaction() const
   {
      if (interactionPtr_ == 0) {
         ExternalPotential* externalPtr = &(system().externalPotential());
         ExternalPotentialImpl< Interaction >* implPtr = 0;
         implPtr = dynamic_cast< ExternalPotentialImpl< Interaction >* >(externalPtr);
         if (implPtr == 0) {
            UTIL_THROW("Failed dynamic cast of ExternalPotential");
         }
         interactionPtr_ = &implPtr->interaction();
      }
      return *interactionPtr_;
   }

   /*
   * Set the external parameter for this McSystem.
   */
   template < class Interaction >
   void McExternalPerturbation<Interaction>::setParameter()
   {  
     interaction().setExternalParameter(parameter_[0]);
   }

   /* 
   * Get the tempering variable from the parent System.
   */
   template < class Interaction >
   double McExternalPerturbation<Interaction>::parameter(int i) const
   { 
      if (i > nParameters_) {
         UTIL_THROW("perturbation parameter index is out of bounds");
      }
     return interaction.externalParameter();
   }

   /* 
   * Return external energy / (kT*epsilon(0,1))
   */
   template < class Interaction >
   double McExternalPerturbation<Interaction>::derivative(int i) const
   {  
      // Preconditions
      if (i > nParameters_) {
         UTIL_THROW("perturbation parameter index is out of bounds");
      }

      if (fabs(parameter_[i] - parameter(i)) > 1.0E-8) {
         UTIL_THROW("Perturbation parameter is not set correctly");
      }
      
      double energy;
      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      
      energy = 0.0;

      for (int iSpec=0; iSpec < system().simulation().nSpecies(); ++iSpec){
           for (system().begin(iSpec, molIter); molIter.notEnd(); ++molIter){
                for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
	            energy += interaction().energy(atomIter->position(), atomIter->typeId());
                    
                }
           }
      }

      // Multiply by the temperature factor.
      if (system().energyEnsemble().isIsothermal()) {
         energy *= system().energyEnsemble().beta();
      } else {
         UTIL_THROW("Non isothermal ensemble.");
      }

      return energy/parameter_[i]; 
   }

}
#endif  // #ifndef SIMP_EXTERNAL
#endif  // #ifdef  MCMD_PERTURB
