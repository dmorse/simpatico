#ifndef EXTERNAL_POTENTIAL_IMPL_H
#define EXTERNAL_POTENTIAL_IMPL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/external/ExternalPotential.h>  // base class
#include <mcMd/simulation/SubSystem.h>                   // base class
#include <util/global.h>

namespace Util
{
   class Vector;
   class Tensor;
}

namespace McMd
{

   using namespace Util;

   class System;

   /**
   * Template implementation of ExternalPotential.
   *
   * \ingroup External_Module
   */
   template <class Interaction>
   class ExternalPotentialImpl : public ExternalPotential, public SubSystem
   {

   public:

      /** 
      * Constructor.
      */
      ExternalPotentialImpl(System& system);

      /** 
      * Constructor (copied from ExternalPotential)
      */
      ExternalPotentialImpl(ExternalPotentialImpl<Interaction>& other);

      /** 
      * Destructor.
      */
      virtual ~ExternalPotentialImpl();

      /**
      * Read param for external potential.
      * 
      * This method reads the external potential Interaction parameter
      * block. Before calling Interaction::readParam(), it passes
      * simulation().nExternalType() to Interaction::setNAtomType().
      */
      virtual void readParam(std::istream& in);

      /// \name Energy, Force, Stress Interactions
      //@{

      /**
      * Sets external parameter
      *
      * \param externalParameter external parameter of system
      */
      virtual void setExternalParameter(double externalParameter);

      /**
      * Returns external parameter
      *
      * \return external parameter
      */
      virtual double externalParameter() const;

      /**
      * Returns external potential energy of a single particle. 
      *
      * \param position atomic position Vector
      * \param typeId   atom type index
      * \return external potential energy
      */
      virtual double energy(const Vector& position, int typeId) const;

      /**
      * Returns force caused by the external potential.
      *
      * \param position  atom position
      * \param typeId    atom type index
      * \param force     force on the atom (on output)
      */
      virtual void getForce(const Vector& position, int typeId, Vector& force) const;

      #if 0
      /**
      * Return external interaction class name (e.g., "TanhCosineExternal").
      */
      virtual std::string interactionClassName() const;
      #endif

      /**
      * Add external force of an Atom to the total force acting on it.
      */
      void addForces();

      /**
      * Return total external energy of this System.
      */
      double energy() const;

      /**
      * Calculate the external energy for one Atom.
      */
      double atomEnergy(const Atom& atom) const;

      //@}

      /**
      * Return external interaction by reference.
      */
      Interaction& interaction();

      /**
      * Return external interaction by const reference.
      */
      const Interaction& interaction() const;

   private:
  
      Interaction* interactionPtr_;

      bool isCopy_;
 
      template <typename T>
      void computeStressImpl(T& stress) const;

   };

}

#include <mcMd/simulation/System.h> 
#include <mcMd/simulation/Simulation.h> 
#include <mcMd/species/Species.h>
#include <util/boundary/Boundary.h> 

#include <util/space/Dimension.h>
#include <util/space/Vector.h>

#include <fstream>

namespace McMd
{

   using namespace Util;

   /* 
   * Default constructor.
   */
   template <class Interaction>
   ExternalPotentialImpl<Interaction>::ExternalPotentialImpl(System& system)
    : ExternalPotential(),
      SubSystem(system),
      interactionPtr_(0),
      isCopy_(false)
   {  interactionPtr_ = new Interaction(); }
 
   /* 
   * Constructor, copy from ExternalPotentialImpl<Interaction>.
   */
   template <class Interaction>
   ExternalPotentialImpl<Interaction>::ExternalPotentialImpl(
                         ExternalPotentialImpl<Interaction>& other)
    : ExternalPotential(),
      SubSystem(other.system()),
      interactionPtr_(&other.interaction()),
      isCopy_(true)
   {}
 
   /* 
   * Destructor. 
   */
   template <class Interaction>
   ExternalPotentialImpl<Interaction>::~ExternalPotentialImpl() 
   {}

   /* 
   * Read parameters from file.
   */
   template <class Interaction>
   void ExternalPotentialImpl<Interaction>::readParam(std::istream &in) 
   {
      // Read only if not a copy.  Do not indent interaction block.
      if (!isCopy_) {
         readBegin(in, "ExternalPotential");
         interaction().setNAtomType(simulation().nAtomType());
         interaction().setBoundary(system().boundary());
         bool nextIndent = false;
         readParamComposite(in, interaction(), nextIndent);
         readEnd(in);
      }
   }
  
   /* 
   * Set external parameter.
   */
   template <class Interaction>
   void ExternalPotentialImpl<Interaction>::setExternalParameter(double externalParameter) 
   { interaction().setExternalParameter(externalParameter); }
 
   /* 
   * Returns external parameter.
   */
   template <class Interaction>
   double ExternalPotentialImpl<Interaction>::externalParameter() const 
   { return interaction().externalParameter(); }
  
   /*
   * Return external energy of an atom.
   */
   template <class Interaction>
   double 
   ExternalPotentialImpl<Interaction>::energy(const Vector& position, int typeId) 
      const
   { return interaction().energy(position, typeId); }

   /*
   * Return external force on an atom.
   */
   template <class Interaction>
   void ExternalPotentialImpl<Interaction>::getForce(const Vector& position, 
                                               int typeId, Vector& force) const
   { interaction().getForce(position, typeId, force); }

   /* 
   * Return total external potential energy.
   */
   template <class Interaction>
   double ExternalPotentialImpl<Interaction>::energy() const
   {
      System::ConstMoleculeIterator molIter;
      Molecule::ConstAtomIterator atomIter;
      double energy = 0.0;
      for (int iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
         for (begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               energy += interaction().energy(atomIter->position(),
                                                    atomIter->typeId());
            }
         }
      }
      return energy;
   }

   /* 
   * Add external forces to total.
   */
   template <class Interaction>
   void ExternalPotentialImpl<Interaction>::addForces()
   {
      Vector  force; 
      System::MoleculeIterator molIter;
      Molecule::AtomIterator   atomIter;
      for (int iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
         for (begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               interaction().getForce(atomIter->position(), 
                                            atomIter->typeId(), force);
               atomIter->force() += force;
            }
         }
      }
   }

   template <class Interaction>
   double ExternalPotentialImpl<Interaction>::atomEnergy(const Atom &atom) const
   {
      double  energy;

      energy = 0.0;
      energy += interaction().energy(atom.position(), atom.typeId());
      return energy;
   }

   template <class Interaction>
   inline Interaction& ExternalPotentialImpl<Interaction>::interaction()
   { return *interactionPtr_; }

   template <class Interaction>
   inline const Interaction& ExternalPotentialImpl<Interaction>::interaction() const
   { return *interactionPtr_; }

   #if 0
   /*
   * Return external interaction class name.
   */
   template <class Interaction>
   std::string ExternalPotentialImpl<Interaction>::interactionClassName() const
   {  return interaction().className(); }
   #endif

}
#endif
