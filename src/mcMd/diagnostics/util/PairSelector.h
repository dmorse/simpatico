#ifndef MCMD_PAIR_SELECTOR_H
#define MCMD_PAIR_SELECTOR_H

#include <util/global.h>
#include <iostream>

namespace McMd
{

   class Atom;

   /**
   * Selection rule for pairs of Atoms.
   *
   * A PairSelector object defines a selection rule for pairs of atoms
   * in other diagnostics that calculate, e.g., RDFs or nonbonded pair 
   * energies. The rule can require that the pair be intramolecular, 
   * intermolecular, or either, and can require that the typeId for 
   * each atom either match a particular value, or that it be left 
   * unconstrained.
   *
   * A PairSelector has a pair type and two atom type indices.
   *
   * The pair type is a PairType enum, which can be INTRA, INTER, or ALL.
   *
   * The atom types ids are integers. A non-negative integer indicates an
   * atom type id. A negative integer indicates that all types are accepted.
   *
   * The bool match() method returns true if a pair satisfies the rule, or 
   * false otherwise.
   */
   class PairSelector
   {

   public:

      /**
      * Type of atom pair, based on identity of parent molecules.
      *
      * Values:
      *  - INTRA : intra-molecular pair
      *  - INTER : inter-molecular pair
      *  - ALL   : accept intra- and inter-molecular pairs
      */
      enum PairType{INTRA, INTER, ALL};

      /**
      * Constructor.
      */
      PairSelector();

      /**
      * Set policy to avoid double counting (true) or to not avoid (false).
      *
      * \param avoidDoubleCounting Policy: true to avoid, false otherwise.
      */
      void setAvoidDoubleCounting(bool avoidDoubleCounting);

      /**
      * Return true if pair of atoms matches the selector policy.
      *
      * \param atom1 first atom in pair.
      * \param atom2 second atom in pair.
      */
      bool match(const Atom& atom1, const Atom& atom2) const;

      /**
      * Return value of pair type.
      */
      PairType pairType() const;

      /**
      * Return value of type Id for atom 1 (-1 means accept all types).
      */
      int atom1TypeId() const;

      /**
      * Return value of type Id for atom 2 (-1 means accept all types).
      */
      int atom2TypeId() const;

      #ifdef UTIL_MPI

      /**
      * Commit associated MPI DataType.
      */
      static void commitMpiType();

      #endif

   private:

      /// Pair type
      PairType pairType_;

      /// Atom type index for atom1.
      int      atom1TypeId_;

      /// Atom type index for atom1.
      int      atom2TypeId_;

      /// Exclude pairs that would otherwise be counted twice in a double loop?
      bool     avoidDoubleCounting_;

   //friends:

      friend std::istream& operator>>(std::istream& in, PairSelector& selector);
      friend std::ostream& operator<<(std::ostream& out, const PairSelector& selector);

   }; 

   /*
   * Return value of pair type.
   */
   inline PairSelector::PairType PairSelector::pairType() const
   {  return pairType_; }

   /*
   * Return value of type Id for atom 1 (-1 means accept all types).
   */
   inline int PairSelector::atom1TypeId() const
   {  return atom1TypeId_; }

   /*
   * Return value of type Id for atom 2 (-1 means accept all types).
   */
   inline int PairSelector::atom2TypeId() const
   {  return atom2TypeId_; }

   // Iostream operator declarations

   /**
   * istream extractor (>>) for a PairSelector object.
   *
   * Format:    pairType  atom1TypeId  atom2TypeId
   *
   * \param in        input stream
   * \param selector  PairSelector to be read from stream
   * \return modified input stream
   */
   std::istream& operator>>(std::istream& in, PairSelector& selector);

   /**
   * ostream inserter (<<) for a PairSelector object.
   *
   * Format:    pairType  atom1TypeId  atom2TypeId
   *
   * \param  out      output stream
   * \param  selector PairSelector to be written to stream
   * \return modified output stream
   */
   std::ostream& operator<<(std::ostream& out, const PairSelector& selector);

   /**
   * istream extractor (>>) for a PairSelector::PairType enum.
   *
   * \param in        input stream
   * \param type      PairType to be read from stream
   * \return modified input stream
   */
   std::istream& operator>>(std::istream& in, PairSelector::PairType& type);

   /**
   * ostream inserter (<<) for a PairSelector::PairType enum.
   *
   * \param  out   output stream
   * \param  type  PairSelector::PairType to be written to stream
   * \return modified output stream
   */
   std::ostream& operator<<(std::ostream& out, const PairSelector::PairType& type);

} 

#ifdef UTIL_MPI
#include <util/mpi/MpiTraits.h>

namespace Util{

   /**
   * Explicit specialization MpiTraits<PairSelector>.
   */
   template <>
   class MpiTraits< McMd::PairSelector>
   {  
   public:  
      static MPI::Datatype type;     ///< MPI Datatype
      static bool hasType;           ///< Is the MPI type initialized?
   };

   /**
   * Explicit specialization MpiTraits<PairSelector::PairType>.
   */
   template <>
   class MpiTraits< McMd::PairSelector::PairType>
   {  
   public:  
      static MPI::Datatype type;     ///< MPI Datatype
      static bool hasType;           ///< Is the MPI type initialized?
   };

}
#endif

#endif
