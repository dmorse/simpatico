#ifndef TOOLS_TYPE_MAP_CPP
#define TOOLS_TYPE_MAP_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TypeMap.h"
#include <util/param/Label.h>

namespace Tools 
{

   using namespace Util;

   /**
   * Constructor
   */
   TypeMap::TypeMap()
   {}

   /**
   * Destructor
   */
   TypeMap::~TypeMap()
   {}

   void TypeMap::insert(int id, const std::string& name)
   {
      std::map<std::string, int>::const_iterator id_iter = ids_.find(name);
      if (id_iter != ids_.end()) {
         UTIL_THROW("Name already present");
      }
      std::map<int, std::string>::const_iterator name_iter = names_.find(id);
      if (name_iter != names_.end()) {
         UTIL_THROW("Id already present");
      }
      ids_.insert(std::pair<std::string, int>(name, id));
      names_.insert(std::pair<int, std::string>(id, name));
   }

   void TypeMap::read(std::istream& in)
   {
      int nType = 0;
      int id;
      std::string name;
      in >> Label("nType") >> nType;
      in >> Label("types");
      for (int i = 0; i < nType; ++i) {
         in >> id >> name;
         insert(id, name);
      }
   }

}
#endif
