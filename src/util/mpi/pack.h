#ifndef PACK_H
#define PACK_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace Util
{

   /*
   * Pack an object of type T into send buffer.
   */
   template <typename T>
   char* pack(char* current, char* end, T* data, int count)
   {
      if (current + sizeof(data) > end) {
         UTIL_THROW("Attempted write past end of send buffer");
      }
      T* ptr = (T*)current;
      for (int i = 0; i < count; ++i) {
         *ptr = data[i];
         ++ptr;
      }
      current = (char *)ptr;
      return current;
   }

   /*
   * Unpack an object of type T from recvBuffer.
   */
   template <typename T>
   char* unpack(char* current, char* end, T* data, int count)
   {
      if (current + sizeof(data) > end) {
         UTIL_THROW("Attempted read past end of recv buffer");
      }
      T* ptr = (T *)current;
      for (int i = 0; i < count; ++i) {
         data[i] = *ptr;
         ++ptr;
      }
      current = (char *)ptr;
      return current;
   }

}
#endif

