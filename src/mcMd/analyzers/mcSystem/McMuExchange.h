#ifndef MCMD_MC_MU_EXCHANGE_H
#define MCMD_MC_MU_EXCHANGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/base/SystemAnalyzer.h>  // base class template
#include <mcMd/mcSimulation/McSystem.h>     // base template parameter
#include <util/accumulators/Average.h>      // member
#include <util/containers/DArray.h>         // member template
#include <util/containers/DSArray.h>        // member template
#include <util/archives/Serializable.h>     // typedef

#include <cstdio>

namespace McMd
{

   using namespace Util;

   class Atom;

   /**
   * Exchange chemical potential for semigrand ensemble.
   *
   * This class calculates the excess free energy cost of
   * transforming a randomly chosen molecule of a specified
   * species into another species by a hypothetical "alchemical" 
   * transformation of the atom types of some or all atoms
   * in the molecule. 
   *
   * See \ref mcMd_analyzer_McMuExchange_page "here" for the
   * parameter file format and any other user documentation.
   *
   * \ingroup McMd_Analyzer_Mc_Module
   */
   class McMuExchange : public SystemAnalyzer<McSystem>
   {

   public:

      /**   
      * Constructor.
      */
      McMuExchange(McSystem& system);

      /**
      * Read parameters and initialize.
      */
      virtual void readParameters(std::istream& in);

      /** 
      * Clear accumulators.
      */
      virtual void setup();
   
      /**
      * Evaluate energy per particle, and add to ensemble. 
      */
      virtual void sample(long iStep);
   
      /**
      * Output results at end of simulation.
      */
      virtual void output();

      /**
      * Save state to binary file archive.
      *
      * \param ar binary saving (output) archive.
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Load state from a binary file archive.
      *
      * \param ar binary loading (input) archive.
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Serialize to/from an archive. 
      *
      * \param ar      saving or loading archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   private:

      /// Get parent McMd::Simulation object.
      Simulation& simulation();

      /// Get Boundary object of parent McSystem.
      Boundary& boundary();

      /// Pointer to Simulation of parent System.
      Simulation  *simulationPtr_;

      /// Pointer to Boundary of parent System.
      Boundary  *boundaryPtr_;

      /// Output file stream
      std::ofstream outputFile_;

      /// Array of accumulators - one per molecule
      DArray<Average> accumulators_;

      /// Array of new atom type ids for all atoms in molecule.
      DArray<int> newTypeIds_;

      /// Array of boolean variables: 1 if atom type changes, 0 if not.
      DArray<int> isAtomFlipped_;

      /// Array of ids for atoms that change type.
      DSArray<int> flipAtomIds_;

      /// Array to hold pointers to neighboring atoms.
      CellList::NeighborArray neighbors_;

      /// Actual number of trial positions for each regrown atom.
      int  speciesId_; 

      /// Number of molecules in this species.
      int  nMolecule_; 

      /// Actual number of trial positions for each regrown atom.
      int  nAtom_; 

      /// Has readParam been called?
      bool isInitialized_;

   };

   // Inline methods

   /*
   * Get Simulation object of parent McSystem.
   */
   inline Simulation& McMuExchange::simulation()
   {  return *simulationPtr_; }

   /*
   * Get Boundary object of parent McSystem.
   */
   inline Boundary& McMuExchange::boundary()
   {  return *boundaryPtr_; }

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void McMuExchange::serialize(Archive& ar, 
                                const unsigned int version)
   {
      ar & interval_;
      ar & outputFileName_;
      ar & speciesId_;
      ar & nAtom_;
      ar & newTypeIds_;
      ar & nMolecule_;
      for (int i = 0; i < nMolecule_; ++i) {
         ar & accumulators_[i];
      }
   }

}
#endif 
