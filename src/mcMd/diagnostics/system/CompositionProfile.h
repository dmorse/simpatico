#ifndef COMPOSITION_PROFILE_H
#define COMPOSITION_PROFILE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
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
   * Evaluates average monomer concentrations in slabs. 
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
      *   - int               nDirections     number of direction vectors
      *   - DArray<Vector>    Directions      direction vectors
      *
      * \param in input parameter stream
      */
      virtual void readParam(std::istream& in);

      /** 
      * Clear accumulators.
      */
      virtual void initialize();
   
      /**
      * Add particle positions to CompositionProfile histogram.
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
      
      /// Array of Miller index IntVectors for wavevectors.
      DArray<IntVector>  intVectors_;

      /// Array of direction vectors.
      DArray<Vector>        waveVectors_;

      /// Number of direction vectors.
      int                      nDirections_;

      /// Number of samples thus far.
      int                      nSample_;

      /// Number of atom types, copied from Simulation::nAtomType().
      int                      nAtomType_;

      /// Has readParam been called?
      bool    isInitialized_;

      /**
      * Update wavevectors.
      */
      void makeWaveVectors();


   };

}
#endif
