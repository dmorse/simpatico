/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "NveVvIntegrator.h"
#include <mcMd/mdSimulation/MdSystem.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/potentials/pair/MdPairPotential.h>
#include <util/boundary/Boundary.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/space/Vector.h>
#include <util/archives/Serializable_includes.h>

//#define USE_ITERATOR 1
namespace McMd
{

   using namespace Util;

   /* 
   * Constructor.   
   */
   NveVvIntegrator::NveVvIntegrator(MdSystem& system)
   : MdIntegrator(system)
   {  setClassName("NveVvIntegrator"); }

   /* 
   * Destructor.   
   */
   NveVvIntegrator::~NveVvIntegrator() 
   {}

   /* 
   * Read parameter and configuration files, initialize system.
   */
   void NveVvIntegrator::readParameters(std::istream &in) 
   {
      read<double>(in, "dt", dt_);
      int nAtomType = simulation().nAtomType();
      prefactors_.allocate(nAtomType);
   }

   /*
   * Load the internal state to an archive.
   */
   void NveVvIntegrator::loadParameters(Serializable::IArchive& ar)
   {  
      loadParameter<double>(ar, "dt", dt_);
      int nAtomType = simulation().nAtomType();
      prefactors_.allocate(nAtomType);
      ar & prefactors_;
   }

   /*
   * Save the internal state to an archive.
   */
   void NveVvIntegrator::save(Serializable::OArchive& ar)
   {  
      ar & dt_;
      ar & prefactors_;
   }

   /* 
   * Initialize constants.
   */
   void NveVvIntegrator::setup() 
   {
      double mass;
      int nAtomType = prefactors_.capacity();
      for (int i = 0; i < nAtomType; ++i) {
         mass = simulation().atomType(i).mass();
         prefactors_[i] = 0.5*dt_/mass;
      }
      system().positionSignal().notify();
      system().velocitySignal().notify();
   }

   /*
   * Verlet MD NVE integrator step
   *
   * This method implements the algorithm:
   *
   *        vm(n)  = v(n) + 0.5*a(n)*dt
   *        x(n+1) = x(n) + vm(n)*dt
   *
   *        calculate acceleration a(n+1)
   *
   *        v(n+1) = vm(n) + 0.5*a(n+1)*dt
   *
   * where x is position, v is velocity, and a is acceleration.
   */
   void NveVvIntegrator::step() 
   {
      Vector dv;
      Vector dr;
      System::MoleculeIterator molIter;
      double prefactor;
      int iSpecies, nSpecies;

      nSpecies = simulation().nSpecies();

      // 1st half velocity Verlet, loop over atoms 
      #if USE_ITERATOR
      Molecule::AtomIterator atomIter;
      for (iSpecies=0; iSpecies < nSpecies; ++iSpecies) {
         system().begin(iSpecies, molIter); 
         for ( ; molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {

               prefactor = prefactors_[atomIter->typeId()];

               dv.multiply(atomIter->force(), prefactor);
               atomIter->velocity() += dv;

               dr.multiply(atomIter->velocity(), dt_);
               atomIter->position() += dr;

            }
         }
      }
      #else 
      Atom* atomPtr;
      int   ia;
      for (iSpecies=0; iSpecies < nSpecies; ++iSpecies) {
         system().begin(iSpecies, molIter); 
         for ( ; molIter.notEnd(); ++molIter) {
            for (ia=0; ia < molIter->nAtom(); ++ia) {
               atomPtr = &molIter->atom(ia);
               prefactor = prefactors_[atomPtr->typeId()];

               dv.multiply(atomPtr->force(), prefactor);
               atomPtr->velocity() += dv;
               dr.multiply(atomPtr->velocity(), dt_);
               
               atomPtr->position() += dr;

            }
         }
      }
      #endif
      system().positionSignal().notify();
      system().velocitySignal().notify();

      #ifndef SIMP_NOPAIR
      if (!system().pairPotential().isPairListCurrent()) {
         system().pairPotential().buildPairList();
      }
      #endif

      system().calculateForces();

      // 2nd half velocity Verlet, loop over atoms
      #if USE_ITERATOR
      for (iSpecies=0; iSpecies < nSpecies; ++iSpecies) {
         system().begin(iSpecies, molIter); 
         for ( ; molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               prefactor = prefactors_[atomIter->typeId()];
               dv.multiply(atomIter->force(), prefactor);
               atomIter->velocity() += dv;
            }
         }
      }
      #else
      for (iSpecies=0; iSpecies < nSpecies; ++iSpecies) {
         system().begin(iSpecies, molIter); 
         for ( ; molIter.notEnd(); ++molIter)
         {
            for (ia=0; ia < molIter->nAtom(); ++ia) {
               atomPtr = &molIter->atom(ia);
               prefactor = prefactors_[atomPtr->typeId()];
               dv.multiply(atomPtr->force(), prefactor);
               atomPtr->velocity() += dv;
            }
         }
      }
      #endif
      system().velocitySignal().notify();

   }

}
#undef USE_ITERATOR
