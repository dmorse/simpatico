/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

#include <util/misc/FlagSet.h>

namespace Util
{

   /*
   * Default constructor.
   */
   FlagSet::FlagSet()
   {}

   /*
   * Constructor.
   */
   FlagSet::FlagSet(std::string allowed)
   {  setAllowed(allowed); }

   /*
   * Set root prefix for all path names.
   */
   void FlagSet::setAllowed(std::string allowed)
   {
      // Copy string
      allowed_ = allowed;

      // Create map, initializing active flag to false for all.
      map_.clear();
      char c;
      int n = allowed_.size();
      for (int i=0; i < n; ++i) {
         c = allowed_[i];
         map_.insert(std::pair<char, bool>(c, false));
      }
   }

   /*
   * Set root prefix for all path names.
   */
   void FlagSet::readOrdered(std::string string)
   {
      // Set all isActive flags to false.
      MapType::iterator iter;
      for (iter = map_.begin(); iter != map_.end(); ++iter) {
         iter->second = false;
      }
 
      // Set flags for characters in string to true.
      int n = allowed_.size();
      if (n > 0) {
         char c, m;
         int j = 0;
         m = string[j];
         for (int i; i < n; ++i) {
            c = string[i];
            while (c != m) {
               ++j;
               if (j == n) {
                  UTIL_THROW("Unknown character");
               }
               m = string[j];
               assert(map_.count(m));
            }
            iter = map_.find(m);
            iter->second = true;
         }
      }
   }

   #if 0
   /*
   */
   void FlagSet::readUnordered(std::string string)
   {
      int n = allowed_.size();
      if (n > 0) {
         int i = 0;
         int j = 0;
         MapType::iterator = iter;
         c = string[i];
         while (i < n) {
            iter = map_.find(c);
            if (iter == MaptType::end) {
               UTIL_THROW("Unknown character");
            }
            
         }
      }
   }
   #endif
}
