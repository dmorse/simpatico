#ifndef SIMP_POINT_H
#define SIMP_POINT_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <simp/species/Species.h>
#include <util/global.h>

namespace Simp
{

   using namespace Util;

   /**
   * A Species in which each Molecule contains only one Atom.
   *
   * \ingroup Simp_Species_Module
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
      * Load species structure from an Archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadSpeciesParam(Serializable::IArchive &ar);

   };
   
   
} 
#endif
