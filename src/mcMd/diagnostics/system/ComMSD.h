#ifndef COM_MSD_H
#define COM_MSD_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, Jian Qin and David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>  // base class template
#include <mcMd/simulation/System.h>             // base class template parameter
#include <util/accumulators/MeanSqDispArray.h>  // member template 
#include <util/space/Vector.h>                  // template parameter
#include <util/space/IntVector.h>               // template parameter
#include <util/containers/DArray.h>             // member template

namespace McMd
{

   using namespace Util;

   class Species;

   /**
   * Molecular center of mass mean squared displacement.
   *
   * \ingroup Diagnostic_Module
   */
   class ComMSD : public SystemDiagnostic<System>
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
      virtual void readParam(std::istream& in);
  
      /** 
      * Determine number of molecules and allocate memory.
      */
      virtual void initialize();
   
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
      if (!isInitialized_) {
         UTIL_THROW("Error: Object not initialized.");
      }

      ar & accumulator_;
      ar & truePositions_;
      ar & oldPositions_;
      ar & shifts_;
      ar & nMolecule_;

      serializeCheck(ar, speciesId_, "speciesId");
      serializeCheck(ar, nAtom_, "nAtom");
      serializeCheck(ar, capacity_, "capacity");

   }

}
#endif
