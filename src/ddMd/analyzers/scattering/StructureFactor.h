#ifndef DDMD_STRUCTURE_FACTOR_H
#define DDMD_STRUCTURE_FACTOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/Analyzer.h>
#include <ddMd/simulation/Simulation.h>
#include <util/containers/DMatrix.h>              // member template
#include <util/containers/DArray.h>               // member template

#include <util/global.h>

#include <iostream>
#include <complex>

namespace DdMd
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
   *     \rho(k) = \sum_{i} v(m, a)\exp( i k \cdot r_i )
   * \f]
   * in which the sum is taken over all particles, \f$a\f$ is the atom type 
   * of particle \f$i\f$, and \f$v(m, a)\f$ is a coefficient for atom type 
   * \f$a\f$ for mode \f$m\f$. The matrix \f$v(m, a)\f$ of mode coefficients 
   * is input the user. 
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
   * \sa \ref ddMd_analyzer_StructureFactor_page "param file format"
   * 
   * \ingroup DdMd_Analyzer_Module
   */
   class StructureFactor : public Analyzer
   {

   public:

      /// Number of analyzer samples in simulation.
      static const int Samples = 100000;
      
      /**	
      * Constructor.
      *
      * \param simulation  reference to parent Simulation object
      */
      StructureFactor(Simulation& simulation);

      /**	
      * Destructor.
      */
      ~StructureFactor();

      /**
      * Read parameters from file.
      *
      * \param in  input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load internal state from an archive.
      *
      * \param ar  input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar  output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);
  
      /** 
      * Clear accumulators.
      */
      virtual void clear();
   
      /**
      * Add particles to StructureFactor accumulators.
      *
      * \param iStep  MD time step counter
      */
      void sample(long iStep);

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
      * First index is wavevector, second index is mode index.
      */
      DMatrix<double> structureFactors_;
      
      /**
      * Fourier modes of concentration.
      *
      * First index is wavevector, second is atom type.
      */
      DMatrix< std::complex<double> >  fourierModes_;

      /**
      * Total fourier modes of concentration.
      *
      * First index is wavevector, second is atom type.
      */
      DMatrix< std::complex<double> >  totalFourierModes_;

      /**
      * Array of Miller index IntVectors for wavevectors.
      */
      DArray<IntVector>  waveIntVectors_;

      /**
      * Array of floating point wave vectors.
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

      /// Has readParam been called?
      bool isInitialized_;

      /// Is this the first step?
      bool isFirstStep_;

      /**
      * Update wavevectors.
      */
      void makeWaveVectors();

   };

}
#endif
