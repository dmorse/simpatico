#ifndef MCMD_VAN_HOVE_H
#define MCMD_VAN_HOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>  // base class template
#include <mcMd/simulation/System.h>             // base class template parameter
#include <util/containers/DArray.h>             // member template
#include <util/accumulators/AutoCorr.h>         // member template parameter
#include <util/space/Vector.h>                  // member template parameter

#include <util/global.h>

#include <iostream>
#include <complex>

namespace McMd
{

   using namespace Util;

   /**
   * Evaluates the van Hove function S(k,t) for one or more wavevector k.
   *
   * The van Hove function S(k,t) is defined here as an expectation value
   * \f[
   *     S(k,t)  = < \psi(k,t) \psi^{*}(k,0) > / V
   * \f]
   * where V is system volume, and \f$\psi(k,t)\f$ is given by a sum
   * \f[
   *     \psi(k,t) = \sum_{i} c_{a} \exp( i k \cdot r_i )
   * \f]
   * over all atoms in the system, where \f$ c_{a} \f$ is a user-specified 
   * coefficient for monomers of type a. 
   *
   * The Van Hove class can calculate S(k,t) for a list of wavevectors.
   * Each wavevector is specified as an IntVector containing integer 
   * Miller indices. These are the coefficients in an expansion of a 
   * reciprocal lattice wavevector in terms of recprocal lattice basis 
   * vectors.
   *
   * In a system with nAtomType = 2, with two monomer types 0 and 1,
   * the parameter file input for calculating the autocorrelation of
   * an order parameter \f$ \psi = \rho_0 - \rho_1 \f$ for 5 wavevectors 
   * along the x axis in an orthorhombic unit cell might look like 
   * this:
   *
   * \code
   * VanHove{
   *    interval                      1000
   *    outputFileName             VanHove
   *    atomTypeCoeffs              1.0000      -1.0000
   *    nBuffer                        100
   *    nWave                            5
   *    waveIntVectors     8      0      0
   *                       9      0      0 
   *                      10      0      0 
   *                      11      0      0 
   *                      12      0      0 
   * }
   * \endcode
   * At the end of a simulation, all of the autocorrelation functions 
   * are output in a file with a suffix *.dat. Each line in this file
   * contains the 3 Miller indices of a wavevector, the absolute
   * magnitude of the wavevector, and a list of nAtomTypeIdPair structure 
   * factor values for the wavevector, one for each atomTypeId pair.
   * 
   * \ingroup McMd_Diagnostic_Module
   */
   class VanHove : public SystemDiagnostic<System>
   {
   public:

      /**	
      * Constructor.
      *
      * \param system reference to parent System object
      */
      VanHove(System &system);

      /**	
      * Destructor.
      */
      ~VanHove();

      /**
      * Read parameters from file.
      *
      * Input format:
      *
      *   - int               interval        sampling interval 
      *   - string            outputFileName  output file base name
      *   - DArray<double>    atomTypeCoeffs  array of type coefficients
      *   - int               nBuffer         number of samples in buffer
      *   - int               nWave           number of wavevectors
      *   - DArray<IntVector> waveIntVectors  IntVector wavevectors
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
      * Save state to an archive.
      *
      * \param ar saving (output) archive.
      */
      virtual void save(Serializable::OArchive& ar);
  
      /**
      * Serialize to/from an archive. 
      * 
      * \param ar      archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /** 
      * Clear accumulators.
      */
      virtual void setup();
   
      /**
      * Add particle pairs to VanHove histogram.
      *
      * \param iStep step counter
      */
      void sample(long iStep);

      /**
      * Output results to predefined output file.
      */
      virtual void output();

   private:

      /// Output file stream.
      std::ofstream outputFile_;

      /// Autocorrelation function accumulators.
      DArray< AutoCorr< std::complex<double>, std::complex<double> >  >
              accumulators_;

      /// Fourier modes. First index is wavevector, second is atom type.
      DArray< std::complex<double> >  fourierModes_;

      /// Array of miller index vectors for wavevectors
      DArray<IntVector>  waveIntVectors_;

      /// Array of wave vectors.
      DArray<Vector>  waveVectors_;

      /// Array of coefficients for atom types.
      DArray<double>  atomTypeCoeffs_;

      /// Number of wavevectors.
      int  nWave_;

      /// Number of samples stored in history buffer.
      int  nBuffer_;

      /// Number of samples thus far.
      int  nSample_;

      /// Number of atom types, copied from Simulation::nAtomType().
      int  nAtomType_;

      /// Has readParam been called?
      bool isInitialized_;

      /// Update wavevectors.
      void makeWaveVectors();

   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void VanHove::serialize(Archive& ar, const unsigned int version)
   {
      Diagnostic::serialize(ar, version);
      ar & nAtomType_;
      ar & atomTypeCoeffs_;
      ar & nBuffer_;
      ar & nWave_;
      ar & waveIntVectors_;
      ar & accumulators_;
      ar & nSample_;
   }

}
#endif
