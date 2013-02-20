#ifndef MCMD_STRUCTURE_FACTOR_P_H
#define MCMD_STRUCTURE_FACTOR_P_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>   // base class template
#include <mcMd/simulation/System.h>              // base class template parameter
#include <util/containers/DMatrix.h>             // member template
#include <util/containers/DArray.h>              // member template
#include <util/containers/Pair.h>                // member template parameter

#include <util/global.h>

#include <iostream>
#include <complex>

namespace McMd
{

   using namespace Util;

   /**
   * StructureFactorP evaluates partial structure factors in Fourier space.
   *
   * A partial structure factor for a wavevector k defined as an expectation 
   * value
   * \f[
   *     S(k)  = < \rho_a(k) \rho_b(-k) / V >
   * \f]
   * where, V is system volume, and \f$\rho_a(k)\f$ is a Fourier component
   * \f[
   *     \rho_a(k) = \sum_{i} \exp( i k \cdot r_i )
   * \f]
   * in which the sum is taken over all particles in a set defined by the
   * particle type index a. 
   *
   * If 0 <= a < nAtomType, then, \f$ \rho_a(k) \f$ is given by a sum over 
   * all atoms of type a. 
   *
   * If a < 0, then \f$ \rho_a(k) \f$ is given by a sum over all atoms, of 
   * all types. By convention, we use a= -1 to indicate this.
   *
   * This class can evaluate multiple types of structures factors,
   * corresponding to different choices for the type pair ab, for
   * multiple wavevectors. The input format thus requires the user
   * to specify one or more atom typeId pairs, and one or more
   * wavevectors. Each wavevector is specified as an IntVector 
   * containing integer Miller indices. These are the coeficients
   * in an expansion of a reciprocal lattice wavevector in terms of
   * recprocal lattice basis vectors for the periodic unit cell.
   *
   * In a system with nAtomType = 2, with two monomer types 0 and 1,
   * the parameter file input for calculating the total, 00 and 01
   * correlation functions for 5 wavevectors along the x axis in an 
   * orthorhombic unit cell might look like this:
   *
   * \code
   * StructureFactorP{
   *    interval                     1000
   *    outputFileName     StructureFactorP
   *    nAtomTypeIdPair                  3
   *    atomTypeIdPairs          -1     -1
   *                              0      0
   *                              0      1
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
   * magnitude of the wavevector, and a list of nAtomTypeIdPair 
   * structure factor values for the wavevector, one for each 
   * atomTypeId pair.
   * 
   * \ingroup McMd_Diagnostic_Module
   */
   class StructureFactorP 
    : public SystemDiagnostic<System>
   {

   public:

      /**	
      * Constructor.
      *
      * \param system reference to parent System object
      */
      StructureFactorP(System &system);

      /**	
      * Destructor.
      */
      ~StructureFactorP();

      /**
      * Read parameters from file.
      *
      * Input format:
      *
      *   - int               interval        sampling interval 
      *   - string            outputFileName  output file base name
      *   - int               nAtomTypeIdPair number of Atom typeId pairs
      *   - DArray<Pair<int>> atomTypeIdPairs integer typeId pairs
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
      * Add particle pairs to StructureFactorP histogram.
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

      /**
      * Structure factor accumulators. 
      * 
      * First index is wavevector, second index is atom typeId pair.
      */
      DMatrix<double>  structureFactors_;

      /*
      * Fourier modes (temporary).
      *
      * First index is wavevector, second is atom type.
      */
      DMatrix< std::complex<double> > fourierModes_;

      /// Array of miller index IntVectors for wavevectors.
      DArray<IntVector>  waveIntVectors_;

      /// Array of wave vectors (temporary)
      DArray<Vector> waveVectors_;

      /// Array of atom type indices (-1 indicates a sum of all types)
      DArray< Pair<int> >  atomTypeIdPairs_;

      /// Number of wavevectors.
      int  nWave_;

      /// Number of selected atom type pairs.
      int  nAtomTypeIdPair_;

      /// Number of samples thus far.
      int  nSample_;

      /// Number of atom types, copied from Simulation::nAtomType().
      int  nAtomType_;

      /// Has readParam been called?
      bool  isInitialized_;

      /**
      * Update wavevectors.
      */
      void makeWaveVectors();

   };

   /**
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void StructureFactorP::serialize(Archive& ar, const unsigned int version)
   {
      Diagnostic::serialize(ar, version);
      ar & nAtomType_;
      ar & nAtomTypeIdPair_;
      ar & atomTypeIdPairs_;
      ar & nWave_;
      ar & waveIntVectors_;

      ar & structureFactors_;
      ar & nSample_;
   }

}
#endif
