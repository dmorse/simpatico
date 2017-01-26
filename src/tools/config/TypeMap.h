#ifndef TOOLS_TYPE_MAP_H
#define TOOLS_TYPE_MAP_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

#include <string>
#include <map>

namespace Tools 
{

   using namespace Util;

   /**
   * Map between type names and type ids.
   *
   * \ingroup Tools_Storage_Module
   */
   class TypeMap 
   {

   public:

      /**
      * Constructor
      */
      TypeMap();

      /**
      * Destructor
      */
      ~TypeMap();

      /**
      * Add typeId / type name pair to the container.
      *
      * \param id  integer type id
      * \param name  type names string
      */
      void insert(int id, const std::string& name);

      /**
      * Read pairs from file.
      *
      * \param in input stream.
      */
      void read(std::istream& in);

      /**
      * Get type id associated with a type name.
      */
      int id(const std::string& name) const;

      /**
      * Get type name associated with an integer id.
      */
      const std::string& name(int id) const;

      /**
      * Number of a pairs in map.
      */
      int size() const;

   private:

      /// Map of string name keys and integer id values.
      std::map<std::string, int> ids_;
 
      /// Map of integer id keys and string name values.    
      std::map<int, std::string> names_;

   };

   inline
   int TypeMap::id(const std::string& name) const
   {
      std::map<std::string, int>::const_iterator iter = ids_.find(name);
      if (iter == ids_.end()) {
         UTIL_THROW("Type name not found");
      }
      return iter->second;   
   }

   inline
   const std::string& TypeMap::name(int id) const
   {
      std::map<int, std::string>::const_iterator iter = names_.find(id);
      if (iter == names_.end()) {
         UTIL_THROW("Type id not found");
      }
      return iter->second;   
   }

   inline 
   int TypeMap::size() const
   {  return names_.size(); }

}
#endif
