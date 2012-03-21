#ifndef MCMD_POINT_H
#define MCMD_POINT_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/species/Species.h>
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /**
   * A Species in which each Molecule contains only one Atom.
   *
   * \ingroup Species_Module
   */
   class Point : public Species
   {
   
   public:
   
      /// Default constructor.
      Point();
   
      /// Destructor.
      virtual ~Point()
      {}
   
   protected:
   
      /**
      * Atom type id for all molecules of this species.
      */
      int type_;
   
      /**
      * Read atom type.
      *
      * \param in input stream
      */
      virtual void readSpeciesParam(std::istream &in);
   
      /**
      * Return the atom type id.
      */
      virtual int getAtomTypeId(Molecule &molecule, int index);
   
      /**
      * Return same bond type for any bond in any chain.
      *
      * \param moleculeId  index of molecule within species
      * \param molBondId   index of bond within the molecule
      */
      //virtual int getBondTypeId(Molecule &molecule, int index);
   
   };
   
   
} 
#endif
