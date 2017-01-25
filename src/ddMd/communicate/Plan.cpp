/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/space/Dimension.h>
#include "Plan.h"

namespace DdMd
{

   unsigned int Plan::GMask[3][2] = { {0x0001, 0x0002}, {0x0004, 0x0008}, {0x0010, 0x0020} };
   unsigned int Plan::EMask[3][2] = { {0x0100, 0x0200}, {0x0400, 0x0800}, {0x1000, 0x2000} };

   using namespace Util;

   /* 
   * Input a Plan from an istream, as a pair of 0s and ones
   */
   std::istream& operator>>(std::istream& in, Plan &plan)
   {
      std::string string;
      int i, j, k;

      in >> string;
      k = 0;
      for (i=0; i < Dimension; ++i) {
         for (j=0; j < 2; ++j) {
            if (string[k] == '1') {
               plan.setExchange(i, j);
            } else if (string[k] == '0') {
               plan.clearExchange(i, j);
            }
            ++k;
         }
      }

      in >> string;
      k = 0;
      for (int i=0; i < Dimension; ++i) {
         for (int j=0; j < 2; ++j) {
            if (string[k] == '1') {
               plan.setGhost(i, j);
            } else if (string[k] == '0') {
               plan.clearGhost(i, j);
            }
         }
         ++k;
      }
      return in;
   }
   
   /* 
   * Output a Plan to an ostream, without line breaks.
   */
   std::ostream& operator<<(std::ostream& out, const Plan &plan) 
   {
      int i, j;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < 2; ++j) {
            if (plan.exchange(i, j)) {
               out << '1';
            } else {
               out << '0';
            }
         }
      }
      out << "  ";
      for (i=0; i < Dimension; ++i) {
         for (j=0; j < 2; ++j) {
            if (plan.ghost(i, j)) {
               out << '1';
            } else {
               out << '0';
            }
         }
      }
      return out;
   }

}
