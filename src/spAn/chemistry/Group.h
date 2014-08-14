#ifndef SPAN_GROUP_H
#define SPAN_GROUP_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace SpAn
{

   /**
   * A group of covalently interacting atoms.
   *
   * \ingroup SpAn_Chemistry_Module
   */
   template <int N>
   class Group
   {
   public: 
   
      /**   
      * Array of integer ids of atoms in this group.
      */
      int atomIds[N];
  
      /** 
      * Integer index for the type of group.
      */
      int typeId;
  
      /** 
      * Global id for this group.
      */
      int id;

   };


   // Associated function declarations

   /**
   * istream extractor (>>) for a Group.
   *
   * \param in        input stream
   * \param group  Group to be read from stream
   * \return modified input stream
   */
   template <int N>
   std::istream& operator>>(std::istream& in, Group<N> &group)
   {
      in >> group.id;
      in >> group.typeId;
      for (int i=0; i < N; ++i) {
         in >> group.atomIds[i];
      }
      return in;
   }
   
   /**
   * ostream inserter (<<) for a Group.
   *
   * Format on one line with no line break:
   *
   * \param  out   output stream
   * \param  group Group to be written to stream
   * \return modified output stream
   */
   template <int N>
   std::ostream& operator << (std::ostream& out, const Group<N> &group) 
   {
      out.width(10);
      out << group.id;
      out.width(10);
      out << group.typeId;
      for (int i = 0; i < N; ++i) {
         out.width(10);
         out << group.atomIds[i];
      }
      return out;
   }

} 
#endif
