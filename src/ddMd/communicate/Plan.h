#ifndef DDMD_PLAN_H
#define DDMD_PLAN_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>

namespace DdMd
{

   /**
   * Communication plan.
   * 
   * A plan contains intructions for a communication pattern for 
   * exchanging atoms, ghosts, or Groups.
   *
   * A Plan has an exchange flag and storage flag for each of the
   * directions in the Plimpton communication scheme. Each direction
   * is indexed by a Cartesian index i=0,..,Dim-1 and j=0,1, where
   * Dim is the dimensionality of space (i.e., normally Dim=3). 
   *
   * When used as a plan for an atom, each exchange flag should be
   * set true in directions for which an atom should be communicated 
   * to exchange ownership, and each ghost flag should be set true 
   * for directions in which the atom must be sent as a ghost in 
   * order to provide position information to other processors.  
   *
   * When used as a plan for a Group, exchange flags are used to
   * denote directions in which the Group must be sent, and ghost
   * flags to denote directions in which the atoms of the Group
   * must be sent as ghosts when the Group is divided among two
   * or more processors. 
   *
   * Implementation: These 12 flags are stored in different bits 
   * of a single unsigned int that can also be accessed or set 
   * directly.
   *
   * \ingroup DdMd_Communicate_Module
   */
   class Plan
   {

   public:

      /**
      * Constructor.
      *
      * All flags cleared upon construction.
      */
      Plan() :
       flags_(0)
      {}

      /**
      * Set exchange flag for direction i, j (set true).
      *
      * \param i Cartesian axis index i=0,1,2=x,y,z
      * \param j binary direction index j=0, 1
      */
      void setExchange(int i, int j)
      {  flags_ |= EMask[i][j]; }
 
      /**
      * Clear exchange flag for direction i, j (set false).
      *
      * \param i Cartesian axis index i=0,1,2=x,y,z
      * \param j binary direction index j=0 (up) 1 (down)
      */
      void clearExchange(int i, int j)
      {  flags_ &= (~EMask[i][j]); }
 
      /**
      * Set ghost flag for direction i, j (set true).
      *
      * \param i Cartesian axis index i=0,1,2=x,y,z
      * \param j binary direction index j=0 (up) 1 (down)
      */
      void setGhost(int i, int j)
      {  flags_ |= GMask[i][j]; }
 
      /**
      * Clear ghost flag for direction i, j (set false).
      *
      * \param i Cartesian axis index i=0,1,2=x,y,z
      * \param j binary direction index j=0 (up) 1 (down)
      */
      void clearGhost(int i, int j)
      {  flags_ &= (~GMask[i][j]); }

      /**
      * Set all flags (contains all bits).
      */
      void setFlags(unsigned int flags)
      {  flags_ = flags; }

      /**
      * Clear all flags (set all to false, set flags_ = 0).
      */ 
      void clearFlags()
      {  flags_ = 0; }
 
      /**
      * Get bool exchange flag for direction i, j.
      *
      * \param i Cartesian axis index i=0,1,2=x,y,z
      * \param j binary direction index j=0 (up) 1 (down)
      */
      bool exchange(int i, int j) const
      {  return bool(flags_ & EMask[i][j]); }
 
      /**
      * Get ghost flag for direction i, j.
      *
      * \param i Cartesian axis index i=0,1,2=x,y,z
      * \param j binary direction index j=0 (up) 1 (down)
      */
      bool ghost(int i, int j) const
      {  return bool(flags_ & GMask[i][j]); }
 
      /**
      * Return raw flags unsigned int.
      */
      unsigned int flags() const
      {  return flags_; }

   private:

      unsigned int flags_;

      /// Matrix of bit masks for exchange flags.
      static unsigned int EMask[3][2];

      /// Matrix of bit masks for ghost flags.
      static unsigned int GMask[3][2];
  
   //friends:

      friend std::istream& operator >> (std::istream& in, Plan &plan);
      friend std::ostream& operator << (std::ostream& out, const Plan &plan);
 
   };

   /**
   * istream extractor (>>) for a Plan.
   *
   * \param in    input stream
   * \param plan  Plan to be read from stream
   * \return modified input stream
   */
   std::istream& operator>>(std::istream& in, Plan &plan);

   /**
   * ostream inserter (<<) for a Plan.
   *
   * Format, one one line with no line break:
   *
   * \param out  output stream
   * \param plan Plan to be written to stream
   * \return modified output stream
   */
   std::ostream& operator<<(std::ostream& out, const Plan &plan);

}
#endif
