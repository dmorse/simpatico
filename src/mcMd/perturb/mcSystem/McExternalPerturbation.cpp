#ifdef  MCMD_PERTURB
#ifdef MCMD_EXTERNAL
#ifndef MC_EXTERNAL_PERTURBATION_CPP
#define MC_EXTERNAL_PERTURBATION_CPP


/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "McExternalPerturbation.h"
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/chemistry/Atom.h>
#include <mcMd/ensembles/EnergyEnsemble.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   McExternalPerturbation::McExternalPerturbation(McSystem& system)
    : LinearPerturbation<McSystem>(system)
   {  }

   /*
   * Destructor.
   */
   McExternalPerturbation::~McExternalPerturbation()
   { }

   /*
   * Read epsilon(0,1) from file
   */
   void McExternalPerturbation::readParam(std::istream& in)
   {  
      Perturbation::readParameters(in); 
      nParameters_ = Perturbation::getNParameters();
   }

   /*
   * Set the parameter epsilon(0,1) for this McSystem.
   */
   void McExternalPerturbation::setParameter()
   {  
     system().externalPotential().setExternalParameter(parameter_[0]);
   }

   /* 
   * Get the tempering variable from the parent System.
   */
   double McExternalPerturbation::parameter(int i) const
   { 
      if (i > nParameters_) {
         UTIL_THROW("perturbation parameter index is out of bounds");
      }
     return system().externalPotential().externalParameter();
   }

   /* 
   * Return pair energy for unlike pairs / (kT*epsilon(0,1))
   */
   double McExternalPerturbation::derivative(int i) const
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
           for (system().begin(iSpec, molIter); !molIter.atEnd(); ++molIter){
                for (molIter->begin(atomIter); !atomIter.atEnd(); ++atomIter) {
	            energy += system().externalPotential().energy(atomIter->position(), atomIter->typeId());
                    
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
#endif  
#endif  // #ifndef MCMD_EXTERNAL
#endif  // #ifdef  MCMD_PERTURB
