#ifndef DDMD_SP_GROUP_H
#define DDMD_SP_GROUP_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace DdMd
{

   /**
   * A group of covalently interacting atoms.
   *
   * \ingroup DdMd_Sp_Chemistry_Module
   */
   template <int N>
   class SpGroup
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
   * istream extractor (>>) for a SpGroup.
   *
   * \param in        input stream
   * \param group  SpGroup to be read from stream
   * \return modified input stream
   */
   template <int N>
   std::istream& operator>>(std::istream& in, SpGroup<N> &group)
   {
      in >> group.id;
      in >> group.typeId;
      for (int i=0; i < N; ++i) {
         in >> group.atomIds[i];
      }
      return in;
   }
   
   /**
   * ostream inserter (<<) for a SpGroup.
   *
   * Format on one line with no line break:
   *
   * \param  out   output stream
   * \param  group SpGroup to be written to stream
   * \return modified output stream
   */
   template <int N>
   std::ostream& operator << (std::ostream& out, const SpGroup<N> &group) 
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
