#ifndef MCMD_LINK_POTENTIAL_IMPL_H
#define MCMD_LINK_POTENTIAL_IMPL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/bond/BondPotential.h>  // base class
#include <mcMd/simulation/SystemInterface.h>           // base class
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
   * an instance of a bond interaction class, but uses the LinkMaster of
   * the parent System to identify linked atoms.
   *
   * \ingroup McMd_Link_Module
   */
   template <class Interaction>
   class LinkPotentialImpl : public BondPotential, public SystemInterface
   {

   public:

      /** 
      * Constructor.
      */
      LinkPotentialImpl(System& system);

      /** 
      * Copy constructor.
      */
      LinkPotentialImpl(LinkPotentialImpl<Interaction>& other);

      /** 
      * Destructor.
      */
      virtual ~LinkPotentialImpl();

      /**
      * Read potential energy.
      * 
      * This method reads the bond potential Interaction parameter
      * block. Before calling Interaction::readParameters(), it passes
      * simulation().nBondType() to Interaction::setNBondType().
      */
      virtual void readParameters(std::istream& in);

      /// \name Bond interaction interface.
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
      * Modify a parameter, identified by a string.
      *
      * \param name   parameter name
      * \param type   bond type index 
      * \param value  new value of parameter
      */
      void set(std::string name, int type, double value)
      {   interactionPtr_->set(name, type, value); }

      /**
      * Get a parameter value, identified by a string.
      *
      * \param name   parameter name
      * \param type bond type index 1
      */
      double get(std::string name, int type) const
      {   return interactionPtr_->get(name, type); }

      /**
      * Return pair interaction class name (e.g., "HarmonicBond").
      */
      virtual std::string interactionClassName() const;

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
      * Compute and increment link forces for all atoms.
      */
      virtual void addForces();

      /**
      * Calculate and store pair energy for this System.
      */
      virtual void computeEnergy();

      /**
      * Compute and store the total nonbonded pressure
      */
      virtual void computeStress();

      //@}

      /**
      * Return bond interaction by reference.
      */
      Interaction& interaction();

      /**
      * Return bond interaction by const reference.
      */
      const Interaction& interaction() const;

   private:
  
      // Pointer to a link/bond interaction.
      Interaction*  interactionPtr_;

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
#include <mcMd/links/Link.h>

#include <simp/species/Species.h>

#include <util/boundary/Boundary.h> 
#include <util/space/Dimension.h>
#include <util/space/Vector.h>
#include <util/space/Tensor.h>
#include <util/accumulators/setToZero.h>

#include <fstream>

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /* 
   * Default constructor.
   */
   template <class Interaction>
   LinkPotentialImpl<Interaction>::LinkPotentialImpl(System& system)
    : BondPotential(),
      SystemInterface(system),
      interactionPtr_(0),
      linkMasterPtr_(&system.linkMaster()),
      isCopy_(false)
   {
      setClassName("LinkPotential"); 
      interactionPtr_ = new Interaction(); 
   }
 
   /* 
   * Constructor, copy from LinkPotentialImpl<Interaction>.
   */
   template <class Interaction>
   LinkPotentialImpl<Interaction>::LinkPotentialImpl(
                         LinkPotentialImpl<Interaction>& other)
    : BondPotential(),
      SystemInterface(other.system()),
      interactionPtr_(&other.interaction()),
      linkMasterPtr_(&other.system().linkMaster()),
      isCopy_(true)
   { setClassName("LinkPotential"); }
 
   /* 
   * Destructor. 
   */
   template <class Interaction>
   LinkPotentialImpl<Interaction>::~LinkPotentialImpl() 
   {}

   /* 
   * Read parameters from file.
   */
   template <class Interaction>
   void LinkPotentialImpl<Interaction>::readParameters(std::istream &in) 
   {
      // Read only if not a copy.  Do not indent interaction block.
      if (!isCopy_) {
         interaction().setNBondType(simulation().nLinkType());
         bool nextIndent = false;
         addParamComposite(interaction(), nextIndent);
         interaction().readParameters(in);
      }
   }
  
   /*
   * Return bond energy for a single pair.
   */
   template <class Interaction>
   double LinkPotentialImpl<Interaction>::energy(double rsq, int iBondType) 
      const
   { return interaction().energy(rsq, iBondType); }

   /*
   * Return force / separation for a single bonded pair.
   */
   template <class Interaction>
   double LinkPotentialImpl<Interaction>::forceOverR(double rsq, int iBondType)
      const
   { return interaction().forceOverR(rsq, iBondType); }

   /*
   * Return force / separation for a single bonded pair.
   */
   template <class Interaction> double 
   LinkPotentialImpl<Interaction>::
      randomBondLength(Random* random, double beta, int bondTypeId) const
   { return interaction().randomBondLength(random, beta, bondTypeId); }

   /*
   * Return link energy for one Atom. 
   */
   template <class Interaction>
   double LinkPotentialImpl<Interaction>::atomEnergy(const Atom &atom) 
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
         energy += interaction().energy(rsq, linkPtr->typeId());
      }

      return energy;
   }

   /* 
   * Add link forces to all atomic forces.
   */
   template <class Interaction>
   void LinkPotentialImpl<Interaction>::addForces()
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
         force *= interaction().forceOverR(rsq, linkPtr->typeId());
         atom0Ptr->force() += force;
         atom1Ptr->force() -= force;
      }
   }

   /* 
   * Return total link energy.
   */
   template <class Interaction>
   void LinkPotentialImpl<Interaction>::computeEnergy()
   {
      double rsq;
      double energy = 0.0;
      Link* linkPtr;
      int iLink, nLink;
      nLink = linkMasterPtr_->nLink();
      for (iLink = 0; iLink < nLink; ++iLink) {
         linkPtr = &(linkMasterPtr_->link(iLink));
         rsq = boundary().distanceSq(linkPtr->atom0().position(), 
                                     linkPtr->atom1().position());
         energy += interaction().energy(rsq, linkPtr->typeId());
      }
      energy_.set(energy);
   }

   /* 
   * Compute the link pressure or stress
   */
   template <class Interaction>
   template <typename T>
   void LinkPotentialImpl<Interaction>::computeStressImpl(T& stress) const
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
         force *= interaction().forceOverR(rsq, linkPtr->typeId());
         incrementPairStress(force, dr, stress);
      }

      // Normalize by volume 
      stress /= boundary().volume();
      normalizeStress(stress);
   }

   /*
   * Compute all short-range pair contributions to stress.
   */
   template <class Interaction>
   void LinkPotentialImpl<Interaction>::computeStress()
   {
      Tensor stress;
      computeStressImpl(stress);

      // Set value of Setable<double> energy_ 
      stress_.set(stress);
   }

   template <class Interaction>
   inline Interaction& LinkPotentialImpl<Interaction>::interaction()
   { return *interactionPtr_; }

   template <class Interaction>
   inline const Interaction& LinkPotentialImpl<Interaction>::interaction() const
   { return *interactionPtr_; }

   /*
   * Return bond interaction class name.
   */
   template <class Interaction>
   std::string LinkPotentialImpl<Interaction>::interactionClassName() const
   {  return interaction().className(); }

}
#endif
