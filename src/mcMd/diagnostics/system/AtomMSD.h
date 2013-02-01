#ifndef MCMD_ATOM_MSD_H
#define MCMD_ATOM_MSD_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>  // base class template
#include <mcMd/simulation/System.h>             // base class template parameter
#include <util/accumulators/MeanSqDispArray.h>  // member template 
#include <util/space/Vector.h>                   // member template parameter
#include <util/containers/DArray.h>             // member template

namespace McMd
{

   using namespace Util;

   /**
   * Autocorrelation for vector separation of two atoms on a molecule.  
   *
   * \ingroup McMd_Diagnostic_Module
   */
   class AtomMSD : public SystemDiagnostic<System>
   {
   
   public:
  
      /**
      * Constructor.
      *
      * \param system reference to parent System
      */
      AtomMSD(System &system);
  
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
   
      /// Output file stream
      std::ofstream outputFile_;

      /// Statistical accumulator
      MeanSqDispArray<Vector> accumulator_;

      /// Array of position vectors, one per molecule of species.
      DArray<Vector>    truePositions_;
   
      /// Array of position vectors, one per molecule of species.
      DArray<Vector>    oldPositions_;
   
      /// Array of shift vectors, one per molecule of species.
      DArray<IntVector> shifts_;
   
      /// Index of relevant Species.
      int     speciesId_;
   
      /// Local index of atom1
      int     atomId_;
   
      /// Number of molecules in the species (must not change).
      int     nMolecule_;
   
      /// Maximum length of each sequence in AutoCorrArray.
      int     capacity_;
   
      /// Has readParam been called?
      int     isInitialized_;
   
   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void AtomMSD::serialize(Archive& ar, const unsigned int version)
   {
      Diagnostic::serialize(ar, version);
      ar & speciesId_;
      ar & atomId_;
      ar & capacity_;
      ar & accumulator_;
      ar & truePositions_;
      ar & oldPositions_;
      ar & shifts_;
      ar & nMolecule_;
   }

}
#endif
