#ifndef MCMD_MD_EWALD_PAIR_POTENTIAL_IMPL_H
#define MCMD_MD_EWALD_PAIR_POTENTIAL_IMPL_H

#include <mcMd/potentials/pair/MdPairPotential.h>
#include <mcMd/potentials/pair/McPairPotentialImpl.h>
#include <mcMd/mdSimulation/MdSystem.h>
#include <mcMd/chemistry/AtomType.h>
#include <util/containers/Array.h>
#include <util/space/Tensor.h>
#include <util/misc/Setable.h>
#include <util/global.h>

#include <mcMd/potentials/coulomb/EwaldInteraction.h>
#include <mcMd/potentials/coulomb/EwaldRSpaceAccumulator.h>

   /**
   * Implementation of a pair potential for a charged system.
   *
   * This class computes forces and energies for all short ranged
   * pair interactions for a charged system, including both   
   * non-Coulomb (e.g., Lennard-Jones) pair interactions and the
   * short range part of the Coulomb interaction in the Ewald method.
   * The addForces() method adds both types of forces to atom forces,
   * but separate accessors are given for non-Coulombic and short
   * range Coulomb contributions to energy and stress.
   */

namespace Util
{
   class Vector;
}

namespace McMd
{

   using namespace Util;

   class MdSystem;
   class EwaldRSpaceAccumulator;
   class EwaldInteraction;
   class MdEwaldPotential;

   template <class Interaction>
   class MdEwaldPairPotentialImpl : public MdPairPotential
   {
      public:
      /**
      * Constructor.
      */
      MdEwaldPairPotentialImpl(MdSystem& system);

      /**
      * Destructor.
      */
      virtual ~MdEwaldPairPotentialImpl();

      /**
      * Read pair potential interaction and pair list blocks.
      *
      * This method reads the pair potential Interaction parameter and
      * PairList blocks, and initializes an internal PairList. Before
      * calling the Interaction::readParameters method, it passes
      * nAtomType to Interaction::setNAtomType().
      *
      * \param in input parameter stream.
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      /// \name Pair Interaction Interface
      //@{

      /**
      * Return non-coulomb pair energy for a single pair.
      *
      * \param rsq       square distance between atoms in pair
      * \param iAtomType atom type index of 1st atom
      * \param jAtomType atom type index of 2nd atom
      * \return energy of pair
      */
      virtual
      double energy(double rsq, int iAtomType, int jAtomType) const;

      /**
      * Return force / separation for a single pair.
      *
      * \param rsq       square distance between atoms in pair
      * \param iAtomType atom type index of 1st atom
      * \param jAtomType atom type index of 2nd atom
      * \return repulsive force (< 0 if attractive) over distance
      */
      virtual
      double forceOverR(double rsq, int iAtomType, int jAtomType) const;

      /**
      * Return maximum cutoff.
      */
      virtual double maxPairCutoff() const;

      /**
      * Modify a parameter, identified by a string.
      *
      * \param name   parameter name
      * \param i      type index of first atom
      * \param j      type index of first atom
      * \param value  new value of parameter
      */
      void set(std::string name, int i, int j, double value)
      {  pairPtr_->set(name, i, j, value); }

      /**
      * Get a parameter value, identified by a string.
      *
      * \param name   parameter name
      * \param i      type index of first atom
      * \param j      type index of first atom
      */
      double get(std::string name, int i, int j) const
      {  return pairPtr_->get(name, i, j); }

      /**
      * Return pair interaction class name (e.g., "LJPair").
      */
      virtual std::string interactionClassName() const;

      //@}
      /// \name Global force and energy calculators
      //@{


       // Force evaluation, which adds both types of pair force.
      virtual void addForces();

      /** 
      * Thermo computes, which computes nonCoulombic and coulombic
      * parts, but stores them in different accumulator variables.
      * The implementation should also mark both accumulators as set.
      */
      virtual void computeEnergy();
      virtual void computeStress();

      //@}

      // Unset both energy accumulators.
      void unsetEnergy()
      { 
         energy_.unset(); 
         rSpaceAccumulatorPtr_->rSpaceEnergy_.unset(); 
      }

      // Unset both stress accumulators.
      void unsetStress()
      { 
         pairStress_.unset(); 
         rSpaceAccumulatorPtr_->rSpaceStress_.unset(); 
      }

      double rSpaceEnergy() const
      { return rSpaceAccumulatorPtr_->rSpaceEnergy_.value(); }

      // Get non-coulombic pair stress.
      Tensor stress()
      { return pairStress_.value(); }

      // Get non-coulombic pair pressure.
      double pressure();

   protected:

      // Non-Coulombic pair interaction
      Interaction* pairPtr_;

   private:

      const Array<AtomType>* atomTypesPtr_;

      MdCoulombPotential* kspacePtr_;
      MdEwaldPotential* ewaldPtr_; 
 
      // Pointers to EwaldRSpaceAccumulator and EwaldInteraction.
      EwaldRSpaceAccumulator*  rSpaceAccumulatorPtr_;
      EwaldInteraction*         ewaldInteractionPtr_;

      // Non-Coulomb pair accumulators
      Setable<Tensor> pairStress_; 

      // Prefactor for coulomb part.
      double fourpiepsi_;

      bool            isCopy_;
   };
}

#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/mdSimulation/MdSystem.h>
#include <mcMd/simulation/stress.h>
#include <mcMd/neighbor/PairIterator.h>
#include <util/boundary/Boundary.h>
#include <util/math/Constants.h>
#include <util/space/Dimension.h>
#include <util/space/Vector.h>
#include <util/accumulators/setToZero.h>
#include <mcMd/potentials/coulomb/MdCoulombPotential.h>
#include <mcMd/potentials/coulomb/MdEwaldPotential.h>
#include <fstream>

namespace McMd
{

   using namespace Util;

   /*
   * Default constructor.
   */
   template <class Interaction>
   MdEwaldPairPotentialImpl<Interaction>::MdEwaldPairPotentialImpl(MdSystem& system)
    : MdPairPotential(system),
      pairPtr_(0),
      atomTypesPtr_(&system.simulation().atomTypes()),
      isCopy_(false)
   {
         // Get pointer to MdCoulombPotential.
         kspacePtr_ = &system.coulombPotential();
 
         // Dynamic cast to a pointer to MdEwaldPotential.
         ewaldPtr_ = dynamic_cast<MdEwaldPotential*>(kspacePtr_);
 
         rSpaceAccumulatorPtr_ = &ewaldPtr_->rSpaceAccumulator_;
 
         ewaldInteractionPtr_  = &ewaldPtr_->ewaldInteraction_;
 
         pairPtr_ = new Interaction;
         // Pass address of MdEwaldPotential to EwaldPair interaction.
         // Note: Uses implicit cast of MdEwaldPotential to its
         // EwaldParameters base class.
   }

   /*
   * Destructor.
   */
   template <class Interaction>
   MdEwaldPairPotentialImpl<Interaction>::~MdEwaldPairPotentialImpl()
   {
      if (pairPtr_ && !isCopy_) {
         delete pairPtr_;
         pairPtr_ = 0;
      }
   }

   template <class Interaction>
   void MdEwaldPairPotentialImpl<Interaction>::readParameters(std::istream& in)
   {
      // Read pair potential parameters only if not a copy.
      if (!isCopy_) {
         pairPtr_->setNAtomType(simulation().nAtomType());
         bool nextIndent = false;
         addParamComposite(*pairPtr_, nextIndent);
         pairPtr_->readParameters(in);
      }

      readParamComposite(in, pairList_);
 
      // Initialize the PairList using larger cutoff between pair interaction and coulomb 
      // interaction. 
      double cutoff = (ewaldInteractionPtr_->rSpaceCutoff() > pairPtr_->maxPairCutoff()) ?
                       ewaldInteractionPtr_->rSpaceCutoff(): pairPtr_->maxPairCutoff();
      pairList_.initialize(simulation().atomCapacity(), cutoff);

      //Initialize prefactor for coulomb part.
      fourpiepsi_ = 1.0 / (4.0*Constants::Pi*ewaldInteractionPtr_->epsilon()) ; 
   }

   /*
   * Load internal state from an archive.
   */
   template <class Interaction>
   void
   MdEwaldPairPotentialImpl<Interaction>::loadParameters(Serializable::IArchive &ar)
   {
      ar >> isCopy_;
      if (!isCopy_) {
         pairPtr_->setNAtomType(simulation().nAtomType());
         bool nextIndent = false;
         addParamComposite(*pairPtr_, nextIndent);
         pairPtr_->loadParameters(ar);
      }
      loadParamComposite(ar, pairList_);

      //Initialize prefactor for coulomb part.
      fourpiepsi_ = 1.0 / (4.0*Constants::Pi*ewaldInteractionPtr_->epsilon()) ; 
   }

   /*
   * Save internal state to an archive.
   */
   template <class Interaction>
   void MdEwaldPairPotentialImpl<Interaction>::save(Serializable::OArchive &ar)
   {
      ar << isCopy_;
      if (!isCopy_) {
         pairPtr_->save(ar);
      }
      pairList_.save(ar);
   }

   /*
   * Return pair energy for a single pair.
   */
   template <class Interaction> double
   MdEwaldPairPotentialImpl<Interaction>::energy(double rsq, int iAtomType, int jAtomType) const
   {  return pairPtr_->energy(rsq, iAtomType, jAtomType); }

   /*
   * Return force / separation for a single pair.
   */
   template <class Interaction> double
   MdEwaldPairPotentialImpl<Interaction>::forceOverR(double rsq,
                                        int iAtomType, int jAtomType) const
   {
      if (rsq < pairPtr_->cutoffSq(iAtomType, jAtomType)) {
         return pairPtr_->forceOverR(rsq, iAtomType, jAtomType);
      } else {
         return 0.0;
      }
   }

   /*
   * Return maximum cutoff.
   */
   template <class Interaction>
   double MdEwaldPairPotentialImpl<Interaction>::maxPairCutoff() const
   { return pairPtr_->maxPairCutoff(); }

   /*
   * Return pair interaction class name.
   */
   template <class Interaction>
   std::string MdEwaldPairPotentialImpl<Interaction>::interactionClassName() const
   {  return pairPtr_->className(); }

   /*
   * Add nonBonded pair forces to atomic forces.
   */
   template <class Interaction>
   void MdEwaldPairPotentialImpl<Interaction>::addForces()
   {
      // Update PairList if necessary
      if (!isPairListCurrent()) {
         buildPairList();
      }

      PairIterator iter;
      Vector       force;
      double       forceOverR;
      double       rsq;
      Atom        *atom0Ptr;
      Atom        *atom1Ptr;
      int          type0, type1;
      double       qProduct;

      // Loop over nonbonded neighbor pairs
      for (pairList_.begin(iter); iter.notEnd(); ++iter) {
         iter.getPair(atom0Ptr, atom1Ptr);
         rsq = boundary().
               distanceSq(atom0Ptr->position(), atom1Ptr->position(),
                          force);
         type0 = atom0Ptr->typeId();
         type1 = atom1Ptr->typeId();
         qProduct = (*atomTypesPtr_)[type0].charge()*(*atomTypesPtr_)[type1].charge();
         if (rsq < pairPtr_->cutoffSq(type0, type1)) {
            forceOverR =  pairPtr_->forceOverR(rsq, type0, type1);
            forceOverR += ewaldInteractionPtr_->rSpaceForceOverR(rsq, qProduct);
            force *= forceOverR;
            atom0Ptr->force() += force;
            atom1Ptr->force() -= force;
         }
      }

   }

   /*
   * Compute and store pair interaction energy.
   * Does nothing if energy is already set.
   */
   template <class Interaction>
   void MdEwaldPairPotentialImpl<Interaction>::computeEnergy()
   {
      // If the pair energy is already known, do nothing and return.
      //if (energy_.isSet()) return;

      // Update PairList if necessary
      if (!isPairListCurrent()) {
         buildPairList();
      }

      // Loop over pairs
      PairIterator iter;
      Atom *atom0Ptr;
      Atom *atom1Ptr;
      double rsq;
      double noncenergy = 0.0;
      double cenergy    = 0.0;
      double qProduct;

      for (pairList_.begin(iter); iter.notEnd(); ++iter) {
         iter.getPair(atom0Ptr, atom1Ptr);

         qProduct =  (*atomTypesPtr_)[atom0Ptr->typeId()].charge();
         qProduct *= (*atomTypesPtr_)[atom1Ptr->typeId()].charge();

         rsq = boundary().
               distanceSq(atom0Ptr->position(), atom1Ptr->position());
         noncenergy += pairPtr_->energy(rsq, atom0Ptr->typeId(), atom1Ptr->typeId());
         cenergy    += ewaldInteractionPtr_->rSpaceEnergy(rsq, qProduct);
      }

      energy_.set(noncenergy); 
      //fourpiepsi_ is prefactor of rpart coulomb energy.
      rSpaceAccumulatorPtr_->rSpaceEnergy_.set(fourpiepsi_ * 0.5 * cenergy);
   }

   /*
   * Add nonBonded pair forces to atomic forces.
   */
   template <class Interaction>
   void MdEwaldPairPotentialImpl<Interaction>::computeStress()
   {

      // If pair stress is already known, do nothing and return.
      //if (stress_.isSet()) return;

      Tensor stress; Vector dr;
      Vector force;
      double rsq;
      double forceOverR        = 0.0;
      double forceOverRcoulomb = 0.0;
      PairIterator iter;
      Atom* atom1Ptr;
      Atom* atom0Ptr;
      int type0, type1;
      double qProduct;

      // Set all elements of stress tensor to zero.
      setToZero(stress);

      // Loop over nonbonded neighbor pairs
      for (pairList_.begin(iter); iter.notEnd(); ++iter) {
         iter.getPair(atom0Ptr, atom1Ptr);
         rsq = boundary().
               distanceSq(atom0Ptr->position(), atom1Ptr->position(), dr);

         type0 = atom0Ptr->typeId();
         type1 = atom1Ptr->typeId();

         if (rsq < pairPtr_->cutoffSq(type0, type1)) {

            //stress from pair potential.
            force  = dr;
            forceOverR =  pairPtr_->forceOverR(rsq, type0, type1);
            force *= forceOverR;
            incrementPairStress(force, dr, stress);
         }

         if (rsq < ewaldInteractionPtr_->rSpaceCutoff()) {

            //charges.
            qProduct =  (*atomTypesPtr_)[atom0Ptr->typeId()].charge();
            qProduct *= (*atomTypesPtr_)[atom1Ptr->typeId()].charge();

            //stress from r-space coulomb potential.
            force = dr;
            forceOverRcoulomb = ewaldInteractionPtr_->rSpaceForceOverR(rsq, qProduct);
            force *= forceOverRcoulomb;
            incrementPairStress(force, dr, stress);
         }
 
      }

      // Normalize by volume
      stress /= boundary().volume();
      normalizeStress(stress);

      pairStress_.set(stress);
   }

}
#endif 
