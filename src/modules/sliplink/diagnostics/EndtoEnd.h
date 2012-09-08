#ifndef END_TO_END_H
#define END_TO_END_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>  // base class template
#include <mcMd/simulation/System.h>             // class template parameter
#include <util/accumulators/Average.h>          // member
#include <util/containers/DArray.h>             // member template
#include <util/space/Vector.h>                   // member template parameter

#include <cstdio> 
#include <cstring> 

namespace McMd
{

   using namespace Util;

   class Species;

   /**
   * End to end distance of a molecule.
   *
   * The radius of gyration R_g is the root mean-squared separation of 
   * atoms in a molecule from the molecular center of mass.  For a species
   * of molecule with N atoms, with positions \f$ R_0, R_1, ...., R_{N-1} \f$,
   * we define
   * \f[
   *  
   *   R_g^2 = \sum_i \langle | R_i - R_{cm} |^2 \rangle / N
   * 
   * \f]
   * where \f$ R_{cm} = (R_0 + R_1 + ... R_{N-1})/N \f$ is the center of mass, 
   * and where \f$ \langle \cdots \rangle \f$ denotes an ensemble average over 
   * all molecules of a specified species. 
   *
   * This class calculates the radius of gyration by accumulating an
   * average value for the above quantity, and optionally outputs block
   * averages to file at an interval specified by the input parameter
   * nSamplePerBlock. No block averages are output if nSamplePerBlock = 0.
   *
   * \ingroup Diagnostic_Module
   */
   class EndtoEnd : public SystemDiagnostic<System>
   {
   
   public:

      /**
      * Constructor.
      *
      * \param system reference to parent System object
      */
      EndtoEnd(System &system);
   
      /**
      * Read parameters from file, and allocate data array.
      *
      * Input format:
      *
      *   - int    interval        : sampling interval
      *   - string outputFileName  : base name for output file(s)
      *   - int    nSamplePerBlock : interval for output of block averages
      *   - int    speciesId       : integer id for Species of interest
      *
      * No block averages are output if nSamplePerBlock = 0. Otherwise,
      * block averages are output to a file named (outputFileName).dat. 
      *
      * \param in parameter input stream
      */
      virtual void readParameters(std::istream& in);
   
      /** 
      * Clear accumulator.
      */
      virtual void setup();
   
      /**
      * Evaluate squared radii of gyration for all molecules, add to ensemble.
      *
      * \param iStep step counter
      */
      virtual void sample(long iStep);
   
      /**
      * Output results at end of simulation.
      */
      virtual void output();

   private:

      /// Output file stream
      std::ofstream outputFile_;

      /// Average object - statistical accumulator
      Average  accumulator_;

      /// Pointer to array of positions for all beads in a molecule
      DArray<Vector> positions_;
   
      /// Pointer to relevant Species.
      Species *speciesPtr_;
   
      /// Number of samples per block average output.
      int nSamplePerBlock_;

      /// Pointer to relevant Species.
      int     speciesId_;
   
      /// Number of atoms per molecule of this species.
      int     nAtom_;
   
   };

}
#endif
