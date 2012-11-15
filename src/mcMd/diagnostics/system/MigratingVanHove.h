#ifdef UTIL_MPI
#ifdef MCMD_PERTURB
#ifndef MCMD_MIGRATING_VAN_HOVE_H
#define MCMD_MIGRATING_VAN_HOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>  // base class template
#include <mcMd/simulation/System.h>             // base class template parameter
#include <mcMd/mcSimulation/McSystem.h>       
#include <util/containers/DArray.h>             // member template
#include <util/accumulators/AutoCorr.h>         // member template parameter
#include <util/space/Vector.h>                  // member template parameter
#include <util/misc/Observer.h>                   
#include <util/misc/Notifier.h>                  
#include <util/mpi/MpiSendRecv.h>
#include <util/mpi/MpiLogger.h>

#include <mcMd/perturb/ReplicaMove.h>          

#include <util/global.h>

#include <iostream>
#include <complex>

namespace McMd
{

   using namespace Util;

   /**
   * MigratingVanHove class evaluates VanHove function 
   * of a configuration by monitoring its trajectory through 
   * different states in space.
   *
   * MigratingVanHove is derived from Observer class (of Observer
   * design pattern) and is therefore notified of any changes in 
   * state by a Notifier (eg. ReplicaExchange, here). Notification is
   * done by invoking the update(sendRecvPair& partners) method,
   * which sends (and receives) all information encapsulated in 
   * VanHove function of current state to (and from) its partner
   * (subsequent) state.
   * 
   * \ingroup McMd_Diagnostic_Module
   */
   class MigratingVanHove : public SystemDiagnostic<System>,
                            public Observer<sendRecvPair>
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
      virtual void readParameters(std::istream& in);

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
      * Send (and receive) VanHove information to (and from) partners.
      *
      * \param partners pair of indices specifying partner states
      */
      virtual void update(const sendRecvPair &partners);

      /**
      * Output results to predefined output file.
      */
      virtual void output();

      /**
      * Save state to archive.
      *
      * \param ar saving (output) archive.
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive.
      */
      virtual void load(Serializable::IArchive& ar);

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
