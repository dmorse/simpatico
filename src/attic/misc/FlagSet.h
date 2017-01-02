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
   * A set of boolean variables represented by characters.
   *
   * A FlagSet has a string of allowed characters, each of which
   * which represents a boolean variable (i.e., a flag), and a
   * string of actual characters containing the subset of the
   * allowed characters that should be set on (i.e., true).
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
      * This function calls setAllowed(string) internally.
      *
      * \param allowed  string of all allowed characters.
      */ 
      FlagSet(std::string allowed);
  
      /**
      * Set or reset the string of allowed flags.
      *
      * This function sets isActive false for all flags and
      * clears the actual string.
      *
      * \param allowed  string of all allowed characters
      */ 
      void setAllowed(std::string allowed);
  
      /**
      * Set the string of actual flag characters.
      *
      * This function requires that the characters in the actual string
      * appear in the same order as they do in the allowed string, but
      * allows some allowed characters to be absent.
      *
      * An Exception is thrown if actual contains a character that is 
      * not allowed, or if it is not in order.
      *
      * \param actual string containing a subset of allowed characters
      */ 
      void setActualOrdered(std::string actual);
  
      /**
      * Is the flag associated with character c active?
      * 
      * \param c  character to be tested.
      */ 
      bool isActive(char c) const;
  
      /**
      * Return the string of allowed characters.
      */ 
      const std::string& allowed() const;
  
      /**
      * Return the string of character for which flags are set.
      */ 
      const std::string& actual() const;
  
   private:

      /// Type of map to allow look up of isActive flags by character.
      typedef std::map<char, bool> MapType;

      /// String containing all allowed characters.
      std::string allowed_;

      /// String containing characters for which flags are set on.
      std::string actual_;

      /// Map containing isActive flags for all characters.
      MapType map_;
 
   };

   // Inline function

   /*
   * Is this flag active?
   */
   inline 
   bool FlagSet::isActive(char c) const
   {
      MapType::const_iterator iter = map_.find(c);
      if (iter == map_.end()) {
         UTIL_THROW("Unknown character");
      } 
      return iter->second;
   }

   /*
   * Return the string of allowed characters.
   */ 
   inline 
   const std::string& FlagSet::allowed() const
   {  return allowed_; }
  
   /*
   * Return the string of actual character flags.
   */ 
   inline 
   const std::string& FlagSet::actual() const
   {  return actual_; }

}
#endif
