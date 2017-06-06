#ifndef MCMD_COM_MSD_H
#define MCMD_COM_MSD_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>  // base class template
#include <mcMd/simulation/System.h>             // base class template parameter
#include <util/accumulators/MeanSqDispArray.h>  // member template 
#include <util/space/Vector.h>                  // template parameter
#include <util/space/IntVector.h>               // template parameter
#include <util/containers/DArray.h>             // member template

namespace Simp {
   class Species;
}

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /**
   * Molecular center of mass mean squared displacement.
   *
   * \sa \ref mcMd_analyzer_ComMSD_page "parameter file format"
   *
   * \ingroup McMd_Analyzer_McMd_Module
   */
   class ComMSD : public SystemAnalyzer<System>
   {
   
   public:
  
      /**
      * Constructor.
      *
      * \param system reference to parent System
      */
      ComMSD(System &system);
  
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
      * Determine number of molecules and allocate memory.
      */
      virtual void setup();
   
      /** 
      * Evaluate center of mass displacement of all chains and add to ensemble.
      *
      * \param iStep counter for number of steps
      */
      virtual void sample(long iStep);
   
      /**
      * Output results to file after simulation is completed.
      */
      virtual void output();

   private:

      /**
      * Compute the center of mass position, given a starting vector.
      *
      * \param R0  position vector of the first atom.
      */
      void moleculeCom(int iMol, Vector &R0);

      /// Output file stream
      std::ofstream outputFile_;

      /// Statistical accumulator
      MeanSqDispArray<Vector> accumulator_;

      /// Array of position vectors, one per molecule of species.
      DArray<Vector>    truePositions_;
   
      /// Array of reference position vectors, one per molecule of species.
      DArray<Vector>    oldPositions_;
   
      /// Array of shift vectors, one per molecule of species.
      DArray<IntVector> shifts_;
   
      /// Index of relevant Species.
      int     speciesId_;
   
      /// Number of molecules in the species (must not change).
      int     nMolecule_;

      /// Number of atoms in per molecule (must not change).
      int     nAtom_;
  
      /// Maximum length of each sequence in AutoCorrArray.
      int     capacity_;

      /// Has readParam been called?
      bool    isInitialized_;

   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void ComMSD::serialize(Archive& ar, const unsigned int version)
   {
      Analyzer::serialize(ar, version);
      ar & speciesId_;
      ar & capacity_;
      ar & nAtom_;
      ar & accumulator_;
      ar & truePositions_;
      ar & oldPositions_;
      ar & shifts_;
      ar & nMolecule_;
   }

}
#endif
