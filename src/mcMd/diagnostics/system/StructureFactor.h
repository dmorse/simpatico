#ifndef MCMD_STRUCTURE_FACTOR_H
#define MCMD_STRUCTURE_FACTOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>    // base class template
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
   * in which the sum is taken over all particles, $a$ is the atom type 
   * of particle $i$, and $v(m, a)$ is a coefficient for atom type $a$
   * for mode $m$. The matrix $v(m, a)$ of mode coefficients is input
   * by the user. 
   *
   * This class can evaluate one or modes for multiple wavevectors.
   * The input format thus requires the user to specify one or more 
   * mode vectors, and one or more wavevectors. Each wavevector is 
   * specified as an IntVector containing the integer Miller indices
   * of a reciprocal lattice vector. These are the coeficients in
   * an expansion of a reciprocal lattice wavevector as a sum of
   * recprocal lattice basis vectors for the periodic unit cell.
   *
   * In a system with nAtomType = 2, with two monomer types 0 and 1,
   * the parameter file input for calculating the density and
   * composition modes for 5 wavevectors along the x axis in an 
   * orthorhombic unit cell might look like this:
   *
   * \code
   * StructureFactor{
   *    interval                      1000
   *    outputFileName    StructureFactor
   *    nMode                            1
   *    modes                     1      1
   *                              0     -1
   *    nWave                            5
   *    waveIntVectors     8      0      0
   *                       9      0      0 
   *                      10      0      0 
   *                      11      0      0 
   *                      12      0      0 
   * }
   * \endcode
   * At the end of a simulation, all of the structure factors are
   * output in a file with a suffix *.dat. Each line in this file
   * contains the 3 Miller indices of a wavevector, the absolute
   * magnitude of the wavevector, and a list of structure factor
   * values for the wavevector, one for each mode.
   * 
   * \ingroup McMd_Diagnostic_Module
   */
   class StructureFactor : public SystemDiagnostic<System>
   {

   public:

      /**
      * Number of diagnostic samples in simulation.
      */
      static const int Samples = 100000;
      
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
      * Serialize to/from an archive. 
      *
      * \param ar      saving or loading archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

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

      /**
      * Array of vector of maximum structure factor values. 
      */
      DArray< std::vector<double> > maximumValue_;

      /**
      * Array of vector of Miller index IntVector with maximum S(q).
      */
      DArray< std::vector<IntVector> > maximumWaveIntVector_;

      /**
      * Array of vector of magnitudes of waveVector with maximum S(q).
      */
      DArray< std::vector<double> > maximumQ_;

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

   private:

      /// Has readParam been called?
      bool isInitialized_;

   };

   /**
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void StructureFactor::serialize(Archive& ar, const unsigned int version)
   {
      Diagnostic::serialize(ar, version);
      ar & nAtomType_;
      ar & nMode_;
      ar & modes_;
      ar & nWave_;
      ar & waveIntVectors_;

      ar & structureFactors_;
      ar & nSample_;

      #if 0
      ar & maximumValue_;
      ar & maximumWaveIntVector_;
      ar & maximumQ_;
      #endif
   }

}
#endif
