#ifndef MCMD_RDF_H
#define MCMD_RDF_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>      // base class template
#include <mcMd/simulation/System.h>                 // base class template parameter
#include <mcMd/diagnostics/util/PairSelector.h>     // member
#include <util/accumulators/RadialDistribution.h>   // member
#include <util/containers/DArray.h>                 // member template

#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /**
   * RDF evaluates the atomic radial distribution function.
   *
   * This class evaluates a radial distribution function in real space,
   * by making a histogram of particle pairs. The algorithm is expensive:
   * It simply executes a double loop over particles, at a cost of order
   * N^2. 
   * 
   * Different types of RDF may be calculated by setting a PairSelector
   * too specify which types of particles pairs should be accepted: The
   * user may specify atomic types for the two particles or accept all 
   * types, and may specify whether to accept only intermolecular,
   * intramolecular or all pairs. Only one histogram (or RDF) can be 
   * calculated.
   *
   * \ingroup McMd_Diagnostic_Module
   */
   class RDF : public SystemDiagnostic<System>
   {
   
   public:

      /**	
      * Constructor.
      *
      * \param system reference to parent System object
      */
      RDF(System &system);

      /**
      * Read parameters from file.  
      *
      * \param in input parameter stream
      */
      virtual void readParam(std::istream& in);
  
      /** 
      * Setup before a simulation (clear accumulator).
      */
      virtual void setup();

      /** 
      * Add particle pairs to RDF histogram.
      *
      * \param iStep step counter
      */
      virtual void sample(long iStep);

      /** 
      * Output results to output file.
      */
      virtual void output();

      /**
      * Save state to binary file archive.
      *
      * \param ar binary saving (output) archive.
      */
      virtual void save(Serializable::OArchiveType& ar);

      /**
      * Load state from a binary file archive.
      *
      * \param ar binary loading (input) archive.
      */
      virtual void load(Serializable::IArchiveType& ar);

      /**
      * Serialize to/from an archive. 
      *
      * \param ar      saving or loading archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   private:

      // Output file stream
      std::ofstream outputFile_;

      // RadialDistribution statistical accumulator
      RadialDistribution  accumulator_;

      /// Sum of snapshot values of number of atoms for each atom type.
      DArray<double> typeNumbers_;

      /// Rule specifying which atom pairs to accept
      PairSelector   selector_;

      /// Sum of snapshot values of concentration for each atom type.
      double normSum_;

      /// Number of atom types, copied from Simulation::nAtomType().
      int nAtomType_;

      /// Has readParam been called?
      bool    isInitialized_;

   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void RDF::serialize(Archive& ar, const unsigned int version)
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: object not initialized");
      }

      // Members that are not set by readParam
      ar & accumulator_;
      ar & typeNumbers_;
      ar & normSum_;

      serializeCheck(ar, nAtomType_, "nAtomType");
   }

}
#endif
