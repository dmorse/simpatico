#ifndef UTIL_FLAG_SET_H
#define UTIL_FLAG_SET_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <string>
#include <map>

namespace Util
{

   /**
   * Represents a set of boolean flags represented by characters.
   *
   * \ingroup Misc_Module
   */
   class FlagSet 
   {

   public:
  
      /**
      * Default constructor.
      */ 
      FlagSet();
  
      /**
      * Constructor.
      *
      * \param allowed  string of all allowed characters.
      */ 
      FlagSet(std::string allowed);
  
      /**
      * Set or reset the string of allowed flags.
      *
      * This function also clears all isActive flags.
      *
      * \param allowed  string of all allowed characters
      */ 
      void setAllowed(std::string allowed);
  
      /**
      * Parse a flag string.
      *
      * \param string  string containing a subset of allowed charaters
      */ 
      void readOrdered(std::string string);
  
      /**
      * Is the flag associated with character c active?
      */ 
      bool isActive(char c) const;
  
   private:

      /// Type of map to allow look up of isActive flags by character.
      typedef std::map<char, bool> MapType;

      /// String containing all allowed characters.
      std::string allowed_;

      /// Map containing isActive flags for all characters.
      MapType map_;
 
   };

   // Inline function

   /*
   * Is this flag active?
   */
   inline bool FlagSet::isActive(char c) const
   {
      MapType::const_iterator iter = map_.find(c);
      if (iter == map_.end()) {
         UTIL_THROW("Unknown character");
      } 
      return iter->second;
   }

}
#endif
