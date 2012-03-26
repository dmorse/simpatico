#ifdef UTIL_MPI
#ifdef MCMD_PERTURB
#ifndef MCMD_MIGRATING_VAN_HOVE_H
#define MCMD_MIGRATING_VAN_HOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>  // base class template
#include <mcMd/simulation/System.h>             // base class template parameter
#include <mcMd/mcSimulation/McSystem.h>             // base class template parameter
#include <util/containers/DArray.h>             // member template
#include <util/accumulators/AutoCorr.h>         // member template parameter
#include <util/space/Vector.h>                   // member template parameter
#include <util/util/Observer.h>                   // member template parameter
#include <util/util/Notifier.h>                   // member template parameter
#include <util/mpi/MpiSendRecv.h>
#include <util/mpi/MpiLogger.h>

#include <mcMd/perturb/ReplicaMove.h>                   // member template parameter

#include <util/global.h>

#include <iostream>
#include <complex>

namespace McMd
{

   using namespace Util;

   /**
   * Evaluates the van Hove function S(k,t) for several k.
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
   class MigratingVanHove : public SystemDiagnostic<System>,
                            public Observer<int>
   {

   public:

      /**	
      * Constructor.
      *
      * \param system reference to parent System object
      */
      MigratingVanHove(System &system);

      /**	
      * Destructor.
      */
      ~MigratingVanHove();

      /**
      * Read parameters from file.
      *
      * Input format:
      *
      *   - int               interval        sampling interval 
      *   - string            outputFileName  output file base name
      *   - DArray<double>    atomTypeCoeffs  array of type coefficients
      *   - int               nBuffer         number of samples in history 
      *   - int               nWave           number of wavevectors
      *   - DArray<IntVector> waveIntVectors  IntVector wavevectors
      *
      * \param in input parameter stream
      */
      virtual void readParam(std::istream& in);

      /** 
      * Clear accumulators.
      */
      virtual void initialize();
   
      /**
      * Add particle pairs to VanHove histogram.
      *
      * \param iStep step counter
      */
      void sample(long iStep);
      
      virtual void update(const int &partnerId);
      /**
      * Output results to predefined output file.
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

      /// Pointer to communicator for replica move.
      MPI::Intracomm* communicatorPtr_;

      /// Pointer to replica move object.
      ReplicaMove* replicaMovePtr_;
     
      /// Has readParam been called?
      bool isInitialized_;

      /// Update wavevectors.
      void makeWaveVectors();

   };


   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void MigratingVanHove::serialize(Archive& ar, const unsigned int version)
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: Object not initialized.");
      }

      ar & accumulators_;
      ar & fourierModes_;
      ar & nSample_;

      serializeCheck(ar, nAtomType_, "nAtomType");
      serializeCheck(ar, nBuffer_, "nBuffer");
      serializeCheck(ar, nWave_, "nWave");
      for (int i = 0; i < nWave_; ++i) {
         serializeCheck(ar, waveIntVectors_[i], "waveIntVector");
      }
   }

}
#endif // ifndef MIGRATING_VAN_HOVE_H
#endif // ifdef UTIL_PERTURB
#endif // ifdef UTIL_MPI
