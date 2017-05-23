#ifndef MCMD_SPECIES_MUTATOR_H
#define MCMD_SPECIES_MUTATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>    // member template
#include <util/param/Label.h>          // member
#include <mcMd/chemistry/Molecule.h>   // inline function

namespace McMd
{

   using namespace Util;

   /**
   * Mix-in class for mutable subclasses of Species.
   *
   * A mutable species is a class that is derived from both Species and
   * SpeciesMutator.  A mutable species models a species in which each
   * molecule can be in any of an enumerable number of internal states.
   * This may be used to represent, for example, homopolymers that can
   * change type within a semi-grand canonical simulation, or polymers 
   * that can grow or shrink within some extended ensemble simulations. 
   * Changes in internal states may change the types of atoms and bonds,
   * and may temporarily erase, or deactivate, certain atoms and bonds.
   *
   * Each molecule in a mutable species is assigned an integer state id, 
   * to specify its internal state. Each state is assigned a floating point 
   * statistical weight. The statistical weights are unnormalized, i.e., 
   * they need not add to unity. The SpeciesMutator class provides arrays 
   * that store these state variables and statistical weights, and methods
   * to access them.
   *
   * \ingroup Simp_Species_Module
   */
   class SpeciesMutator 
   {
   
   public:
   
      /**
      * Constructor.
      *
      * The constructor sets Species::isSpeciesMutator_ to true. 
      */
      SpeciesMutator(); 

      /**
      * Destructor.
      */
      virtual ~SpeciesMutator();

      /// \name Mutators
      //@{

      /**
      * Read the state id for one molecule from a configuration file stream
      * 
      * \param in input stream
      * \param molecule molecule of interest
      */
      virtual void readMoleculeState(std::istream& in, Molecule& molecule);

      /**
      * Write the state id for one molecule to a configuration file stream
      * 
      * \param out      output stream
      * \param molecule molecule of interest
      */
      virtual 
      void writeMoleculeState(std::ostream& out, const Molecule& molecule) const;

      /**
      * Change the state of a specific molecule.
      *
      * The method should be used to implement MC moves that change the internal
      * state of a molecule, after the previous stateId and bead positions are
      * known.  The method can use knowledge of the previous stateId and previous 
      * bead positions to implement transitions between a specific pairs of states.
      *
      * \param molecule reference to molecule of interest.
      * \param stateId  integer index of new state.
      */
      virtual void setMoleculeState(Molecule& molecule, int stateId) = 0;

      /**
      * Set the statistical weight associated with a specific state.
      *
      * \param stateId  integer index of state
      * \param weight   desired statistical weight 
      */
      void setWeight(int stateId, double weight);

      //@{
      // Accessors

      /**
      * Get the state id for a specific molecule.
      *
      * \param molecule reference to molecule of interest.
      */
      int moleculeStateId(const Molecule& molecule) const;

      /**
      * Get the number of possible molecule states.
      */     
      int nState() const;

      /**
      * Get the number of molecules with a specified state.
      *
      * \param stateId index for specific state.
      */
      int stateOccupancy(int stateId) const;

      /**
      * Get the statistical weight for a specfic molecular state.
      *
      * \param stateId index for specific state.
      */
      double stateWeight(int stateId) const;

      //@}

      /**
      * Null value for a state index.
      */
      static const int NullStateId = -1;

   protected:

      /**
      * Allocate arrays of molecule state ids and statistical weights.
      *
      * \param nMolecule number of molecules allocated for species.
      * \param nState    number of possible internal states.
      */
      void allocateSpeciesMutator(int nMolecule, int nState);

      /**
      * Set the state id of a specific molecule.
      *
      * This method should be called by the public method setMoleculeState.
      * to set the state id for a molecule. It does not modify the molecule.
      *
      * \sa setMoleculeState
      *
      * \param molecule reference to molecule of interest.
      * \param stateId  integer index of new state.
      */
      void setMoleculeStateId(const Molecule& molecule, int stateId);

   private:

      /**
      * Array of statistical weights for internal states, indexed by state id.
      *
      * The capacity of stateWeights_ must equal nState_.
      */
      DArray<double>  stateWeights_;

      /**
      * Histogram of number of molecules in each state.
      *
      * The capacity of stateOccupancies_ must equal nState_.
      */
      DArray<int>  stateOccupancies_;

      /**
      * Array of state ids for individial molecules, indexed by molecule id.
      *
      * The capacity must be set equal to the number of molecules in this 
      * species. If A is a Molecule object with an index moleculeId = A.id(), 
      * then moleculeStateIds_[moleculeId] is the state index for object A.
      */
      DArray<int>  moleculeStateIds_;

      /**
      * Label for stateId in a config file.
      */
      Label  stateIdLabel_;

      /**
      * The number of possible different molecular states.
      */
      int  nState_;

   };

   // Inline method definitions
   
   /*
   * Get the state id for a molecule.
   */
   inline 
   int SpeciesMutator::moleculeStateId(const Molecule& molecule) const
   {  return moleculeStateIds_[molecule.id()]; }

   /*
   * Get the occupancy (number of molecules) in a specific internal state.
   */
   inline 
   int SpeciesMutator::stateOccupancy(int stateId) const
   {  return stateOccupancies_[stateId]; }

   /*
   * Get the number of possible internal states.
   */
   inline 
   int SpeciesMutator::nState() const
   {  return nState_; }

   /*
   * Get the statistical weight for a specific internal state.
   */
   inline 
   double SpeciesMutator::stateWeight(int stateId) const
   {  return stateWeights_[stateId]; }

   /*
   * Set the statistical weight for a specific internal state.
   */
   inline 
   void SpeciesMutator::setWeight(int stateId, double weight)
   {  stateWeights_[stateId] = weight; }

} 
#endif
