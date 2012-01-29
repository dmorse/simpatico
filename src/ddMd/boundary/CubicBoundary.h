#ifndef CUBIC_BOUNDARY_H
#define CUBIC_BOUNDARY_H

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
   * A cubic periodic simulation cell.
   *
   * \ingroup Boundary_Module
   */
   class CubicBoundary : public OrthoBoundaryBase
   {

   public:

      /**
      * Constructor.
      */
      CubicBoundary();

      /**
      * Set unit cell dimension.
      *
      * Also sets all related lengths and volume.
      *
      * \param length unit cell length along x, y, or z.
      */
      void setLength(double length);

      /**
      * Return true if valid, or throw Exception.
      */
      bool isValid();

   // friends:

      /// Unit test
      friend class ::OrthoBoundaryTest;

      /// ostream inserter
      friend std::ostream& operator << (std::ostream& out, const CubicBoundary& boundary);

      /// istream extractor
      friend std::istream& operator >> (std::istream& in, CubicBoundary& boundary);

   };

   /**
   * istream extractor for a CubicBoundary.
   *
   * \param  in       input stream
   * \param  boundary CubicBoundary to be read from stream
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, CubicBoundary& boundary);

   /**
   * ostream inserter for a CubicBoundary.
   *
   * \param  out      output stream
   * \param  boundary CubicBoundary to be written to stream
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, const CubicBoundary& boundary);

} 
#endif
