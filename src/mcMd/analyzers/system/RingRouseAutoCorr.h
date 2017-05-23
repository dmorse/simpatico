#ifndef MCMD_RING_ROUSE_AUTO_CORR_H
#define MCMD_RING_ROUSE_AUTO_CORR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>   // base class template
#include <mcMd/simulation/System.h>          // base class template param
#include <util/accumulators/AutoCorrArray.h> // member template 
#include <util/space/Vector.h>               // member template parameter
#include <util/containers/DArray.h>          // member template

namespace Simp {
   class Species;
}

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /**
   * Autocorrelation for Rouse mode coefficients of a ring molecule.
   *
   * \ingroup McMd_Analyzer_McMd_Module
   */
   class RingRouseAutoCorr : public SystemAnalyzer<System>
   {
   
   public:
  
      /**
      * Constructor.
      *
      * \param system reference to parent System
      */
      RingRouseAutoCorr(System &system);

      /**
      * Destructor.
      */
      virtual ~RingRouseAutoCorr();
  
      /** 
      * Read parameters from file, and allocate data array.
      *
      * \param in input parameter stream
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

      /**
      * Serialize internal state to/from archive.
      *
      * \param ar       archive
      * \param version  id for archive version
      */
      template <class Archive>
      void serialize(Archive &ar, const unsigned int version);

      /** 
      * Allocate memory and initialize Rouse mode eigenvector p.
      */
      virtual void setup();

      /** 
      * Evaluate end-to-end vectors of all chains, add to ensemble.
      *
      * \param iStep counter for number of steps
      */
      virtual void sample(long iStep);

      /**
      * Output results to file after simulation is completed.
      */
      virtual void output();
   
   private:
   
      /// Output file stream.
      std::ofstream outputFile_;

      /// Statistical accumulator.
      AutoCorrArray<Vector, double> accumulator_;

      /// Array of Rouse mode coefficients, one per molecule of species.
      DArray<Vector> data_;
   
      /// Array of coefficients for projecting to the Rouse mode (eigenvector).
      DArray<double> projector_;
   
      /// Pointer to relevant Species.
      Species *speciesPtr_;
   
      /// Index of relevant Species.
      int     speciesId_;
   
      /// Number of molecules in the species (must not change).
      int     nMolecule_;

      /// Number of atoms per molecule (must not change).
      int     nAtom_;
 
      /// Index of Rouse mode.
      int     p_;

      /// Maximum length of each sequence in AutoCorrArray.
      int     capacity_;
   
      /// Has readParam been called?
      bool    isInitialized_;

   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void RingRouseAutoCorr::serialize(Archive& ar, const unsigned int version)
   {
      Analyzer::serialize(ar, version);
      ar & speciesId_;
      ar & p_;
      ar & capacity_;
      ar & nAtom_;
      ar & nMolecule_;
      ar & accumulator_;
      ar & projector_;
   }

}
#endif
