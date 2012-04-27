#ifndef MCMD_MAXIMUM_SF_H
#define MCMD_MAXIMUM_SF_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "StructureFactor.h"
#include <util/crystal/LatticeSystem.h>

namespace McMd
{

   using namespace Util;

   /**
   * MaximumSf outputs the wavevector at which S(q) is maximum 
   * and the S(q) value itself.
   *
   * This class evaluates the maximum S(q) by checking through
   * all wavevectors within a star and all stars within a region 
   * in which allof the h, k, l integer wavevector components 
   * (i.e., Miller indices) has an absolute magnitude less than 
   * or equal to a parameter hMax.
   * 
   * At each interval step, the value of maximum S(q) and 
   * the maximum q vector is output in a file with a suffix *.dat. 
   * Each line in this file contains the 3 Miller indices of 
   * wavevector, the absolute magnitude of the wavevector, and 
   * the maximum structure factor value.
   * 
   * \ingroup McMd_Diagnostic_Module
   */
   class MaximumSf : public StructureFactor
   {

   public:

      /**	
      * Constructor.
      *
      * \param system reference to parent System object
      */
      MaximumSf(System &system);

      /**
      * Read parameters from file.
      *
      * Input format:
      *
      *   - int               interval        sampling interval 
      *   - string            outputFileName  output file base name
      *   - int               nMode           number of modes
      *   - DMatrix<double>   modes           mode vectors
      *   - int               hMax            maximum Miller index
      *
      * \param in input parameter stream
      */
      virtual void readParam(std::istream& in);

      /**
      * Set up before a simulation.
      */
      virtual void setup();

      /**
      * Output maximum structure factor at end of simulation.
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
 
   private:

      /// Array of ids for first wavevector in each star.
      DArray<int>  starIds_;
      
      /// Array of star sizes.
      DArray<int>  starSizes_;
      
      /// Maximum Miller index of wavevectors in grid.
      int   hMax_;
      
      /// Number of stars of symmetry related wavevectors.
      int   nStar_;
     
      /// Lattice system used to create stars.
      LatticeSystem   lattice_;
      
      /// Output file stream
      std::ofstream outputFile_;
      
      /// Has readParam been called?
      bool    isInitialized_;
   };

}
#endif
