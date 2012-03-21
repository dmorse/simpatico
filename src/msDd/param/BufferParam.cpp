#ifndef MSDD_BUFFER_PARAM_CPP
#define MSDD_BUFFER_PARAM_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "BufferParam.h"

namespace MsDd
{
   using namespace Util;

   /*
   * Constructor.
   */
   BufferParam::BufferParam()
    : Util::ParamComposite(),
      atomCapacity_(-1),
      ghostCapacity_(-1)
   {}

   /*
   * Destructor.
   */
   BufferParam::~BufferParam()
   {}

   /*
   * Read capacities, and allocate buffers.
   */
   void BufferParam::readParam(std::istream& in)
   {
      readBegin(in, "Buffer");
      read<int>(in, "atomCapacity",  atomCapacity_);
      read<int>(in, "ghostCapacity", ghostCapacity_);
      readEnd(in);
   }

   /*
   * Maximum number of atoms for which space is available.
   */
   int BufferParam::atomCapacity() const
   {  return atomCapacity_; }

   /*
   * Maximum number of ghost atoms for which space is available.
   */
   int BufferParam::ghostCapacity() const
   {  return ghostCapacity_; }

}
#endif
