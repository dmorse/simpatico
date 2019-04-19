#ifndef MCMD_STRUCTURE_FACTOR_H
#define MCMD_STRUCTURE_FACTOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/base/SystemAnalyzer.h>    // base class template
#include <mcMd/simulation/System.h>               // base class template parameter
#include <util/containers/DMatrix.h>              // member template
#include <util/containers/DArray.h>               // member template

#include <util/global.h>

#include <iostream>
#include <complex>

namespace McMd
{

   using namespace Util;

   /**
   * StructureFactor evaluates structure factors in Fourier space.
   *
   * A structure factor for a wavevector k and mode m defined as an 
   * expectation value
   * \f[
   *     S_{m}(k)  = < \rho(k) \rho(-k) / V >
   * \f]
   * where, V is system volume, and \f$\rho(k)\f$ is a Fourier component
   * \f[
   *     \rho(k) = \sum_{a} \sum_{i} v(m, a)\exp( i k \cdot r_i )
   * \f]
   * in which the sum is taken over all particles, a is the atom type 
   * of particle i, and v(m, a) is a coefficient for atom type a for
   * mode m. The matrix v(m, a) of mode coefficients is input by the
   * user. 
   *
   * This class can evaluate one or modes for multiple wavevectors.
   * The input format thus requires the user to specify one or more 
   * mode vectors, and one or more wavevectors. Each wavevector is 
   * specified as an IntVector containing the integer Miller indices
   * of a reciprocal lattice vector. These are the coeficients in
   * an expansion of a reciprocal lattice wavevector as a sum of
   * recprocal lattice basis vectors for the periodic unit cell.
   *
   * \sa \ref mcMd_analyzer_StructureFactor_page "parameter file format"
   * 
   * \ingroup McMd_Analyzer_McMd_Module
   */
   class StructureFactor : public SystemAnalyzer<System>
   {

   public:

      /**	
      * Constructor.
      *
      * \param system reference to parent System object
      */
      StructureFactor(System &system);

      /**	
      * Destructor.
      */
      ~StructureFactor();

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
      * Clear accumulators.
      */
      virtual void setup();
   
      /**
      * Add particles to StructureFactor accumulators.
      *
      * \param iStep step counter
      */
      virtual void sample(long iStep);

      /**
      * Output results to predefined output file.
      */
      virtual void output();

   protected:

      /**
      * Output file stream.
      */
      std::ofstream outputFile_;

      /**
      * Structure factor accumulators. 
      * 
      * First index is wavevector, second is a mode index.
      */
      DMatrix<double> structureFactors_;
      
      /**
      * Instantaneous Fourier amplitudes (temporary)
      *
      * First index is wavevector, second is mode index.
      */
      DMatrix< std::complex<double> > fourierModes_;

      /**
      * Array of Miller index IntVectors for wavevectors.
      */
      DArray<IntVector>  waveIntVectors_;

      /**
      * Array of floating point wave vectors (temporary).
      */
      DArray<Vector>  waveVectors_;

      /**
      * Array of mode vectors 
      *
      * First index is mode, second is atomType.
      */
      DMatrix<double>  modes_;

      /// Number of wavevectors.
      int  nWave_;

      /// Number of mode vectors
      int  nMode_;

      /// Number of samples thus far.
      int  nSample_;

      /// Number of atom types, copied from Simulation::nAtomType().
      int  nAtomType_;

      /**
      * Update wavevectors.
      */
      void makeWaveVectors();

      /// Is this the first step?
      bool isFirstStep_;

   private:

      /// Has readParam been called?
      bool isInitialized_;
   };

}
#endif
