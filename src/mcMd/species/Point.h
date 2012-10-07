#ifndef MCMD_POINT_H
#define MCMD_POINT_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
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
   * \ingroup McMd_Species_Module
   */
   class Point : public Species
   {
   
   public:
   
      /// Default constructor.
      Point();
   
      /// Destructor.
      virtual ~Point()
      {}
   
      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

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
      * Load species structure from an Archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadSpeciesParam(Serializable::IArchive &ar);

   };
   
   
} 
#endif
