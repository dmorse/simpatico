#ifndef DDMD_PLAN_H
#define DDMD_PLAN_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>

namespace DdMd
{

   /**
   * Plan for communication pattern for exchanging atoms
   * and communicating ghots.
   *
   * \ingroup DdMd_Communicate_Module
   */
   class Plan
   {

   public:

      Plan() :
       flags_(0)
      {}

      void setExchange(int i, int j)
      {  flags_ |= EMask[i][j]; }
 
      void clearExchange(int i, int j)
      {  flags_ &= (~EMask[i][j]); }
 
      void setGhost(int i, int j)
      {  flags_ |= GMask[i][j]; }
 
      void clearGhost(int i, int j)
      {  flags_ &= (~GMask[i][j]); }

      void setFlags(unsigned int flags)
      {  flags_ = flags; }
 
      void clearFlags()
      {  flags_ = 0; }
 
      bool exchange(int i, int j) const
      {  return bool(flags_ & EMask[i][j]); }
 
      bool ghost(int i, int j) const
      {  return bool(flags_ & GMask[i][j]); }
 
      unsigned int flags() const
      {  return flags_; }

      // friends:
      friend std::istream& operator >> (std::istream& in, Plan &plan);
      friend std::ostream& operator << (std::ostream& out, const Plan &plan);
 
   private:

      unsigned int flags_;
 
      static unsigned int GMask[3][2];
      static unsigned int EMask[3][2];

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
   * \param  out  output stream
   * \param  plan Plan to be written to stream
   * \return modified output stream
   */
   std::ostream& operator<<(std::ostream& out, const Plan &plan);

}
#endif
