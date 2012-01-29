#ifdef  MCMD_LINK
#ifndef LINK_POTENTIAL_IMPL_H
#define LINK_POTENTIAL_IMPL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/bond/BondPotential.h>  // base class
#include <mcMd/simulation/SubSystem.h>           // base class
#include <mcMd/links/LinkMaster.h>               // member
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
   * Template implementation of an BondPotential for links.
   *
   * A LinkPotential implements the BondPotential interface, and uses
   * an instance of a bond evaluator class, but uses the LinkMaster of
   * the parent System to identify linked atoms.
   *
   * \ingroup Link_Module
   */
   template <class Evaluator>
   class LinkPotentialImpl : public BondPotential, public SubSystem
   {

   public:

      /** 
      * Constructor.
      */
      LinkPotentialImpl(System& system);

      /** 
      * Copy constructor.
      */
      LinkPotentialImpl(LinkPotentialImpl<Evaluator>& other);

      /** 
      * Destructor.
      */
      virtual ~LinkPotentialImpl();

      /**
      * Read potential energy.
      * 
      * This method reads the bond potential Evaluator parameter
      * block. Before calling Evaluator::readParam(), it passes
      * simulation().nBondType() to Evaluator::setNBondType().
      */
      virtual void readParam(std::istream& in);

      /// \name Bond evaluator interface.
      //@{

      /**
      * Return pair energy for a single pair.
      */
      virtual double energy(double rsq, int bondTypeId) const;

      /**
      * Return force / separation for a single pair.
      */
      virtual double forceOverR(double rsq, int bondTypeId) const;

      /**
      * Return force / separation for a single pair.
      */
      virtual 
      double randomBondLength(Random* random, double beta, int bondTypeId) 
      const;

      /**
      * Return pair evaluator class name (e.g., "HarmonicBond").
      */
      virtual std::string evaluatorClassName() const;

      //@}
      /// \name System energy, force and stress.
      //@{

      /**
      * Calculate the link potential energy for one Atom.
      *
      * \param  atom Atom object of interest
      * \return bond energy of one atom. 
      */
      virtual double atomEnergy(const Atom& atom) const;

      /**
      * Return total link potential energy of the System.
      */
      virtual double energy() const;

      /**
      * Compute and increment link forces for all atoms.
      */
      virtual void addForces();

      /**
      * Compute total nonbonded pressure
      *
      * \param stress (output) pressure.
      */
      virtual void computeStress(double& stress) const;

      /**
      * Compute x, y, z nonbonded pressures.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Vector& stress) const;

      /**
      * Compute stress tensor.
      *
      * \param stress (output) pressures.
      */
      virtual void computeStress(Util::Tensor& stress) const;

      //@}

      /**
      * Return bond evaluator by reference.
      */
      Evaluator& evaluator();

      /**
      * Return bond evaluator by const reference.
      */
      const Evaluator& evaluator() const;

   private:
  
      // Pointer to a link/bond evaluator.
      Evaluator*  evaluatorPtr_;

      // Pointer to LinkMaster of parent System.
      LinkMaster* linkMasterPtr_;

      // What this created with a copy constructor?
      bool isCopy_;

      // Generic implementation of stress calculation
      template <typename T>
      void computeStressImpl(T& stress) const;

   };

}

#include <mcMd/simulation/System.h> 
#include <mcMd/simulation/Simulation.h> 
#include <mcMd/simulation/stress.h>
#include <mcMd/species/Species.h>
#include <mcMd/boundary/Boundary.h> 
#include <mcMd/links/Link.h>

#include <util/space/Dimension.h>
#include <util/space/Vector.h>
#include <util/space/Tensor.h>
#include <util/accumulators/setToZero.h>

#include <fstream>

namespace McMd
{

   using namespace Util;

   /* 
   * Default constructor.
   */
   template <class Evaluator>
   LinkPotentialImpl<Evaluator>::LinkPotentialImpl(System& system)
    : BondPotential(),
      SubSystem(system),
      evaluatorPtr_(0),
      linkMasterPtr_(&system.linkMaster()),
      isCopy_(false)
   { evaluatorPtr_ = new Evaluator(); }
 
   /* 
   * Constructor, copy from LinkPotentialImpl<Evaluator>.
   */
   template <class Evaluator>
   LinkPotentialImpl<Evaluator>::LinkPotentialImpl(
                         LinkPotentialImpl<Evaluator>& other)
    : BondPotential(),
      SubSystem(other.system()),
      evaluatorPtr_(&other.evaluator()),
      linkMasterPtr_(&other.system().linkMaster()),
      isCopy_(true)
   {}
 
   /* 
   * Destructor. 
   */
   template <class Evaluator>
   LinkPotentialImpl<Evaluator>::~LinkPotentialImpl() 
   {}

   /* 
   * Read parameters from file.
   */
   template <class Evaluator>
   void LinkPotentialImpl<Evaluator>::readParam(std::istream &in) 
   {
      // Read only if not a copy.  Do not indent evaluator block.
      if (!isCopy_) {
         readBegin(in, "LinkPotential");
         evaluator().setNBondType(simulation().nLinkType());
         bool nextIndent = false;
         readParamComposite(in, evaluator(), nextIndent);
         readEnd(in);
      }
   }
  
   /*
   * Return bond energy for a single pair.
   */
   template <class Evaluator>
   double LinkPotentialImpl<Evaluator>::energy(double rsq, int iBondType) 
      const
   { return evaluator().energy(rsq, iBondType); }

   /*
   * Return force / separation for a single bonded pair.
   */
   template <class Evaluator>
   double LinkPotentialImpl<Evaluator>::forceOverR(double rsq, int iBondType)
      const
   { return evaluator().forceOverR(rsq, iBondType); }

   /*
   * Return force / separation for a single bonded pair.
   */
   template <class Evaluator> double 
   LinkPotentialImpl<Evaluator>::
      randomBondLength(Random* random, double beta, int bondTypeId) const
   { return evaluator().randomBondLength(random, beta, bondTypeId); }

   /*
   * Return link energy for one Atom. 
   */
   template <class Evaluator>
   double LinkPotentialImpl<Evaluator>::atomEnergy(const Atom &atom) 
   const
   {
      double rsq;
      double energy = 0.0;
      LinkMaster::AtomLinkSet* linkSetPtr;
      Link*                    linkPtr;
      int iLink, nLink;

      linkSetPtr = &linkMasterPtr_->atomLinkSet(atom);
      nLink = linkSetPtr->size();
      for (iLink = 0; iLink < nLink; ++iLink) {
         linkPtr = &((*linkSetPtr)[iLink]);
         rsq = boundary().distanceSq(linkPtr->atom0().position(),
                                     linkPtr->atom1().position());
         energy += evaluator().energy(rsq, linkPtr->typeId());
      }

      return energy;
   }

   /* 
   * Return total link energy.
   */
   template <class Evaluator>
   double LinkPotentialImpl<Evaluator>::energy() const
   {
      double rsq;
      double energy = 0.0;
      Link*  linkPtr;
      int    iLink, nLink;
      nLink = linkMasterPtr_->nLink();
      for (iLink = 0; iLink < nLink; ++iLink) {
         linkPtr = &(linkMasterPtr_->link(iLink));
         rsq = boundary().distanceSq(linkPtr->atom0().position(), 
                                     linkPtr->atom1().position());
         energy += evaluator().energy(rsq, linkPtr->typeId());
      }
      return energy;
   }

   /* 
   * Add link forces to all atomic forces.
   */
   template <class Evaluator>
   void LinkPotentialImpl<Evaluator>::addForces()
   {
      Vector  force;
      double  rsq;
      Link*   linkPtr;
      Atom*   atom0Ptr;
      Atom*   atom1Ptr;
      int     iLink, nLink;
      nLink = linkMasterPtr_->nLink();
      for (iLink = 0; iLink < nLink; ++iLink) {
         linkPtr  = &linkMasterPtr_->link(iLink);
         atom0Ptr = &(linkPtr->atom0());
         atom1Ptr = &(linkPtr->atom1());
         rsq = boundary().
               distanceSq(atom0Ptr->position(), atom1Ptr->position(), force);
         force *= evaluator().forceOverR(rsq, linkPtr->typeId());
         atom0Ptr->force() += force;
         atom1Ptr->force() -= force;
      }
   }

   /* 
   * Compute the link pressure or stress
   */
   template <class Evaluator>
   template <typename T>
   void LinkPotentialImpl<Evaluator>::computeStressImpl(T& stress) const
   {
      Vector  dr;
      Vector  force;
      double  rsq;
      Link*   linkPtr;
      int     iLink, nLink;

      setToZero(stress);

      // Loop over links
      nLink = linkMasterPtr_->nLink();
      for (iLink = 0; iLink < nLink; ++iLink) {
         linkPtr  = &linkMasterPtr_->link(iLink);
         rsq = boundary().distanceSq(linkPtr->atom0().position(), 
                                     linkPtr->atom1().position(), dr);
         force  = dr;
         force *= evaluator().forceOverR(rsq, linkPtr->typeId());
         incrementPairStress(force, dr, stress);
      }

      // Normalize by volume 
      stress /= boundary().volume();
      normalizeStress(stress);
   }

   template <class Evaluator>
   void LinkPotentialImpl<Evaluator>::computeStress(double& stress) const
   {  computeStressImpl(stress); }

   template <class Evaluator>
   void LinkPotentialImpl<Evaluator>::computeStress(Util::Vector& stress) 
        const
   {  computeStressImpl(stress); }

   template <class Evaluator>
   void LinkPotentialImpl<Evaluator>::computeStress(Util::Tensor& stress) 
        const
   {  computeStressImpl(stress); }

   template <class Evaluator>
   inline Evaluator& LinkPotentialImpl<Evaluator>::evaluator()
   { return *evaluatorPtr_; }

   template <class Evaluator>
   inline const Evaluator& LinkPotentialImpl<Evaluator>::evaluator() const
   { return *evaluatorPtr_; }

   /*
   * Return bond potential evaluator class name.
   */
   template <class Evaluator>
   std::string LinkPotentialImpl<Evaluator>::evaluatorClassName() const
   {  return evaluator().className(); }

}
#endif
#endif
