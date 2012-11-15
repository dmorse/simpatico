#ifndef MCMD_BOND_LENGTH_DIST_H
#define MCMD_BOND_LENGTH_DIST_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>    // base class template
#include <mcMd/simulation/System.h>               // base class template parameter
#include <util/accumulators/Distribution.h>
#include <util/containers/DArray.h>

#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /**
   * BondLengthDist evaluates the distribution function of the lengths of the bonds.
   *
   * \ingroup McMd_Diagnostic_Module
   */
   class BondLengthDist : public SystemDiagnostic<System>
   {
   
   public:

      /**	
      * Constructor.
      *
      * \param system reference to parent System object
      */
      BondLengthDist(System &system);

      /**
      * Read parameters from file.  
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive.
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save state to archive.
      *
      * \param ar saving (output) archive.
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Serialize to/from an archive. 
      *
      * \param ar      saving or loading archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /** 
      * Clear accumulator.
      */
      virtual void setup();
   
      /** 
      * Add particle pairs to BondLengthDist histogram.
      *
      * \param iStep step counter
      */
      void sample(long iStep);

      /** 
      * Output results to output file.
      */
      virtual void output();

   private:

      // Output file stream
      std::ofstream outputFile_;

      // Distribution statistical accumulator
      Distribution  accumulator_;

      /// Minimum of range
      double min_;

      /// Maximum of range
      double max_;

      /// Number of bins in range
      double nBin_;

      /// Index of relevant Species.
      int     speciesId_;
  
      /// Has readParam been called?
      bool    isInitialized_;

   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void BondLengthDist::serialize(Archive& ar, const unsigned int version)
   {   
      Diagnostic::serialize(ar, version);
      ar & speciesId_;
      ar & min_;
      ar & max_;
      ar & nBin_;
      ar & accumulator_;
   }

}
#endif
