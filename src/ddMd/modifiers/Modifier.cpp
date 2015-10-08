/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Modifier.h"

namespace DdMd
{

   using namespace Util;

   const Bit Modifier::Flags::Setup = 0;
   const Bit Modifier::Flags::PreIntegrate1 = 1;
   const Bit Modifier::Flags::PostIntegrate1 = 2;
   const Bit Modifier::Flags::PreTransform = 3;
   const Bit Modifier::Flags::PreExchange = 4;
   const Bit Modifier::Flags::PostExchange = 5;
   const Bit Modifier::Flags::PostNeighbor = 6;
   const Bit Modifier::Flags::PreUpdate = 7;
   const Bit Modifier::Flags::PostUpdate = 8;
   const Bit Modifier::Flags::PreForce = 9;
   const Bit Modifier::Flags::PostForce = 10;
   const Bit Modifier::Flags::EndOfStep = 11;
   const Bit Modifier::Flags::Exchange = 12;
   const Bit Modifier::Flags::Update = 13;
   const Bit Modifier::Flags::ReverseUpdate = 14;

   /*
   * This static method exists to guarantee initialization of static 
   * constants that are defined in the same file.  Call it somewhere 
   * in the program to guarantee that the contents of this file will 
   * be linked, rather than optimized away.
   */
   void Modifier::initStatic()
   {  
      // This function can only be called once.
      static int nCall = 0;
      if (nCall == 0) {
         ++nCall;
      }
   }

   /*
   * Default constructor.
   */
   Modifier::Modifier()
    : ParamComposite(),
      flags_(0),
      interval_(1),
      simulationPtr_(0)
   {}

   /*
   * Constructor.
   */
   Modifier::Modifier(Simulation& simulation)
    : ParamComposite(),
      flags_(0),
      interval_(1),
      simulationPtr_(&simulation)
   {}

   /*
   * Destructor.
   */
   Modifier::~Modifier()
   {}

   /*
   * Return true if a flag is set, false otherwise.
   */
   bool Modifier::isSet(Bit flag) const
   {  return flag.isSet(flags_); }

   /*
   * Return unsigned int representation of all bit flags.
   */
   unsigned int Modifier::flags() const
   {  return flags_; }

   // Protected methods

   /*
   * Return true if a flag is set, false otherwise.
   */
   void Modifier::set(Bit flag)
   {  flag.set(flags_); }

   /*
   * Read the interval from parameter file, with error checking.
   */
   void Modifier::readInterval(std::istream &in)
   {

      // Read interval value (inherited from Interval)
      read<long>(in, "interval", interval_);

      // Check that interval has a nonzero, positive value
      if (interval_ == 0) {
         UTIL_THROW("interval_ == 0");
      }
      if (interval_ < 0) {
         UTIL_THROW("interval_ < 0");
      }

   }

   /*
   * Load parameter interval from input archive.
   */
   void Modifier::loadInterval(Serializable::IArchive &ar)
   {
      loadParameter<long>(ar, "interval", interval_);

      // Check that interval has a nonzero, positive value
      if (interval_ == 0) {
         UTIL_THROW("interval_ == 0");
      }
      if (interval_ < 0) {
         UTIL_THROW("interval_ < 0");
      }
   }

   /*
   * Save interval parameter to an archive.
   */
   void Modifier::saveInterval(Serializable::OArchive &ar)
   {  ar << interval_; }

}
