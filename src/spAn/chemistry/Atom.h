#ifndef SPAN_ATOM_H
#define SPAN_ATOM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/space/Vector.h>            // members

namespace SpAn
{

   using namespace Util;

   /**
   * A point particle in an MD simulation.
   *
   * \ingroup SpAn_Chemistry_Module
   */
   struct Atom
   {

   public:

      /// Atom position
      Vector position;    

      /// Atom velocity
      Vector velocity;    

      /// Atom type index
      int typeId;         

      /// Unique global index (tag)
      int id;             

      /// Index for species of parent molecule
      int speciesId;      

      /// Index of molecule with its species
      int moleculeId;     

      /// Index of atom within its molecule
      int atomId;         

      /*
      * The speciesId, moleculeId and atomId indices are not contained
      * in all configuration file formats, and are thus optional.
      */

   };

}
#endif
