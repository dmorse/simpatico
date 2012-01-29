#ifndef TETRAGONAL_BOUNDARY_H
#define TETRAGONAL_BOUNDARY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "OrthoBoundaryBase.h"         // base class
#include <iostream>

class OrthoBoundaryTest;

namespace DdMd
{

   using namespace Util;

   /**
   * An tetragonal periodic simulation cell.
   *
   * \ingroup Boundary_Module
   */
   class TetragonalBoundary : public OrthoBoundaryBase
   {

   public:

      /**
      * Constructor.
      */
      TetragonalBoundary();

      /**
      * Set unit cell dimensions.
      *
      * Also sets all related lengths and volume.
      *
      * \param ab unit cell length along x and y axes
      * \param c  unit cell length along the z axis
      */
      void setLengths(double ab, double c);

      /**
      * Return true if valid, or throw Exception.
      */
      bool isValid();

   // friends:

      /// Unit test
      friend class ::OrthoBoundaryTest;

      /// ostream inserter
      friend std::ostream& operator << (std::ostream& out, const TetragonalBoundary& boundary);

      /// istream extractor
      friend std::istream& operator >> (std::istream& in, TetragonalBoundary& boundary);

   };

   /**
   * istream extractor for a TetragonalBoundary.
   *
   * \param  in       input stream
   * \param  boundary TetragonalBoundary to be read from stream
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, TetragonalBoundary& boundary);

   /**
   * ostream inserter for a TetragonalBoundary.
   *
   * \param  out      output stream
   * \param  boundary TetragonalBoundary to be written to stream
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, const TetragonalBoundary& boundary);

} 
#endif
