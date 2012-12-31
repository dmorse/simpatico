#ifndef DDMD_COMPOSITION_PROFILE_H
#define DDMD_COMPOSITION_PROFILE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/diagnostics/SystemDiagnostic.h>    // base class template
#include <ddMd/simulation/Simulation.h>               // base class template parameter
#include <util/accumulators/Distribution.h>      
#include <util/containers/DArray.h>               // member template

#include <util/global.h>

#include <iostream>

namespace DdMd
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
   class CompositionProfile : public Diagnostic
   {

   public:

      /**	
      * Constructor.
      *
      * \param system reference to parent System object
      */
      CompositionProfile(Simulation &simulation);

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
      virtual void clear();
   
      /**
      * Add particle positions to histogram.
      *
      * \param iStep step counter
      */
      void sample(long iStep);

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

}
#endif
