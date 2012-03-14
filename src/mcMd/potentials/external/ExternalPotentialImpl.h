#ifdef  INTER_EXTERNAL
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
   template <class Evaluator>
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
      ExternalPotentialImpl(ExternalPotentialImpl<Evaluator>& other);

      /** 
      * Destructor.
      */
      virtual ~ExternalPotentialImpl();

      /**
      * Read param for external potential.
      * 
      * This method reads the external potential Evaluator parameter
      * block. Before calling Evaluator::readParam(), it passes
      * simulation().nExternalType() to Evaluator::setNAtomType().
      */
      virtual void readParam(std::istream& in);

      /// \name Energy, Force, Stress Evaluators
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

      /**
      * Return external evaluator class name (e.g., "TanhCosineExternal").
      */
      virtual std::string evaluatorClassName() const;

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
      * Return external evaluator by reference.
      */
      Evaluator& evaluator();

      /**
      * Return external evaluator by const reference.
      */
      const Evaluator& evaluator() const;

   private:
  
      Evaluator* evaluatorPtr_;

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
   template <class Evaluator>
   ExternalPotentialImpl<Evaluator>::ExternalPotentialImpl(System& system)
    : ExternalPotential(),
      SubSystem(system),
      evaluatorPtr_(0),
      isCopy_(false)
   {  evaluatorPtr_ = new Evaluator(); }
 
   /* 
   * Constructor, copy from ExternalPotentialImpl<Evaluator>.
   */
   template <class Evaluator>
   ExternalPotentialImpl<Evaluator>::ExternalPotentialImpl(
                         ExternalPotentialImpl<Evaluator>& other)
    : ExternalPotential(),
      SubSystem(other.system()),
      evaluatorPtr_(&other.evaluator()),
      isCopy_(true)
   {}
 
   /* 
   * Destructor. 
   */
   template <class Evaluator>
   ExternalPotentialImpl<Evaluator>::~ExternalPotentialImpl() 
   {}

   /* 
   * Read parameters from file.
   */
   template <class Evaluator>
   void ExternalPotentialImpl<Evaluator>::readParam(std::istream &in) 
   {
      // Read only if not a copy.  Do not indent evaluator block.
      if (!isCopy_) {
         readBegin(in, "ExternalPotential");
         evaluator().setNAtomType(simulation().nAtomType());
         evaluator().setBoundary(system().boundary());
         bool nextIndent = false;
         readParamComposite(in, evaluator(), nextIndent);
         readEnd(in);
      }
   }
  
   /* 
   * Set external parameter.
   */
   template <class Evaluator>
   void ExternalPotentialImpl<Evaluator>::setExternalParameter(double externalParameter) 
   { evaluator().setExternalParameter(externalParameter); }
 
   /* 
   * Returns external parameter.
   */
   template <class Evaluator>
   double ExternalPotentialImpl<Evaluator>::externalParameter() const 
   { return evaluator().externalParameter(); }
  
   /*
   * Return external energy of an atom.
   */
   template <class Evaluator>
   double 
   ExternalPotentialImpl<Evaluator>::energy(const Vector& position, int typeId) 
      const
   { return evaluator().energy(position, typeId); }

   /*
   * Return external force on an atom.
   */
   template <class Evaluator>
   void ExternalPotentialImpl<Evaluator>::getForce(const Vector& position, 
                                               int typeId, Vector& force) const
   { evaluator().getForce(position, typeId, force); }

   /* 
   * Return total external potential energy.
   */
   template <class Evaluator>
   double ExternalPotentialImpl<Evaluator>::energy() const
   {
      System::ConstMoleculeIterator molIter;
      Molecule::ConstAtomIterator atomIter;
      double energy = 0.0;
      for (int iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
         for (begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               energy += evaluator().energy(atomIter->position(),
                                                    atomIter->typeId());
            }
         }
      }
      return energy;
   }

   /* 
   * Add external forces to total.
   */
   template <class Evaluator>
   void ExternalPotentialImpl<Evaluator>::addForces()
   {
      Vector  force; 
      System::MoleculeIterator molIter;
      Molecule::AtomIterator   atomIter;
      for (int iSpec=0; iSpec < simulation().nSpecies(); ++iSpec) {
         for (begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               evaluator().getForce(atomIter->position(), 
                                            atomIter->typeId(), force);
               atomIter->force() += force;
            }
         }
      }
   }

   template <class Evaluator>
   double ExternalPotentialImpl<Evaluator>::atomEnergy(const Atom &atom) const
   {
      double  energy;

      energy = 0.0;
      energy += evaluator().energy(atom.position(), atom.typeId());
      return energy;
   }

   template <class Evaluator>
   inline Evaluator& ExternalPotentialImpl<Evaluator>::evaluator()
   { return *evaluatorPtr_; }

   template <class Evaluator>
   inline const Evaluator& ExternalPotentialImpl<Evaluator>::evaluator() const
   { return *evaluatorPtr_; }

   /*
   * Return external potential evaluator class name.
   */
   template <class Evaluator>
   std::string ExternalPotentialImpl<Evaluator>::evaluatorClassName() const
   {  return evaluator().className(); }

}
#endif
#endif
