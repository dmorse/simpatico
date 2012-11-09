#ifndef MCMD_COMPOSITION_PROFILE_H
#define MCMD_COMPOSITION_PROFILE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>    // base class template
#include <mcMd/simulation/System.h>               // base class template parameter
#include <util/accumulators/Distribution.h>      
#include <util/containers/DArray.h>               // member template

#include <util/global.h>

#include <iostream>

namespace McMd
{

   using namespace Util;
   /**
   * CompositionProfile evaluates the distribution of monomer 
   * positions along several user-specified directions. A direction 
   * vector is specified as an IntVector containing integer 
   * Miller indices. 

   * The dot product of monomer position vector and unit 
   * direction vector is added to distribution function of 
   * particular monomer type and direction vector. 
   */
   class CompositionProfile : public SystemDiagnostic<System>
   {

   public:

      /**	
      * Constructor.
      *
      * \param system reference to parent System object
      */
      CompositionProfile(System &system);

      /**	
      * Destructor.
      */
      ~CompositionProfile();

      /**
      * Read parameters from file.
      *
      * Input format:
      *
      *   - int               interval        sampling interval 
      *   - string            outputFileName  output file base name
      *   - int               nDirection      number of directions
      *   - DArray<IntVector> intVectors      IntVector directions
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /** 
      * Clear accumulators.
      */
      virtual void setup();
   
      /**
      * Add particle positions to histogram.
      *
      * \param iStep step counter
      */
      void sample(long iStep);

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
      * \param version  index for archive version
      */
      template <class Archive>
      void serialize(Archive &ar, const unsigned int version);


      /**
      * Output results to predefined output file.
      */
      virtual void output();

   protected:

      /// Output file stream.
      std::ofstream outputFile_;

      /// Distribution statistical accumulators.
      DArray<Distribution> accumulators_;
      
      /// Array of Miller index vectors for directions.
      DArray<IntVector> intVectors_;

      /// Array of direction vectors.
      DArray<Vector> waveVectors_;

      /// Number of directions.
      int  nDirection_;

      /// Number of samples thus far.
      int  nSample_;

      /// Number of atom types, copied from Simulation::nAtomType().
      int  nAtomType_;

      /// Has readParam been called?
      bool isInitialized_;

      /**
      * Update wavevectors.
      */
      void makeWaveVectors();


   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void CompositionProfile::serialize(Archive& ar, const unsigned int version)
   {
      Diagnostic::serialize(ar, version);
      ar & nDirection_;
      ar & intVectors_;
      ar & waveVectors_;
      ar & accumulators_;
      ar & nSample_;
      ar & nAtomType_;
   }

}
#endif
