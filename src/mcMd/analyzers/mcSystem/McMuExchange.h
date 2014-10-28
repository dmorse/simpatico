#ifndef MCMD_MC_MU_EXCHANGE_H
#define MCMD_MC_MU_EXCHANGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>  // base class template
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
   * \ingroup McMd_Analyzer_Module
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
      * Clear accumulator.
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
      virtual void load(Serializable::IArchive& ar);

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

      /// Average object - statistical accumulator
      Average  accumulator_;

      /// Array of new atom type ids for all atoms in molecule.
      DArray<int> newTypeIds_;

      /// Array of old atom type ids for all atoms in molecule.
      DArray<int> oldTypeIds_;

      /// Array of ids for atoms that change type.
      DSArray<int> flipAtomIds_;

      /// Actual number of trial positions for each regrown atom.
      int  speciesId_; 

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
      if (!isInitialized_) {
         UTIL_THROW("Error: Object not initialized.");
      }

      ar & accumulator_;
   }

}
#endif 
