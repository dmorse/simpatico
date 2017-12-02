/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Atom.h"
#include <ddMd/chemistry/Mask.h>
#include <ddMd/communicate/Plan.h>
#ifdef UTIL_MPI
#include <ddMd/communicate/Buffer.h>
#endif

namespace DdMd
{

   using namespace Util;

   /*
   * Initialize hasAtomContext_ flag to false. This disables storage and
   * communication of AtomContext information by default.
   */
   bool Atom::hasAtomContext_ = false;

   /*
   * Enable (true) or disable (false) use of AtomContext data.
   */
   void Atom::setHasAtomContext(bool hasAtomContext)
   {  hasAtomContext_ = hasAtomContext; }

   /*
   * Constructor (private, used only by AtomArray).
   */
   Atom::Atom() :
     position_(0.0),
     typeId_(-1),
     localId_(0),
     force_(0.0),
     arrayPtr_(0)
     #ifdef UTIL_32BIT
     , pad_(0)
     #endif
   {}

   /*
   * Assignment (public).
   */
   Atom& Atom::operator= (const Atom& other)
   {
      position_ = other.position_;
      typeId_ = other.typeId_;
      force_ = other.force_;
      setIsGhost(other.isGhost());
      velocity() = other.velocity();
      setId(other.id());
      plan() = other.plan();
      groups() = other.groups();
      if (hasAtomContext_) {
         context() = other.context();
      }
      mask() = other.mask();
      return *this;
   }

   /*
   * Reset integer members to null values.
   */
   void Atom::clear()
   {
      typeId_ = -1;
      setIsGhost(false);
      setId(-1);
      plan().clearFlags();
      groups() = 0;
      if (hasAtomContext_) {
         context().clear();
      }
      mask().clear();
   }

   #ifdef UTIL_MPI
   // Atom Exchange

   /*
   * Pack a local Atom for exchange of ownership.
   */
   void Atom::packAtom(Buffer& buffer)
   {
      buffer.pack<int>(id());
      buffer.pack<int>(typeId());
      buffer.pack<Vector>(position());
      buffer.pack<Vector>(velocity());
      buffer.pack<unsigned int>(plan().flags());
      buffer.pack<unsigned int>(groups());
      if (hasAtomContext_) {
         buffer.pack<AtomContext>(context());
      }

      // Pack Mask
      Mask& m = mask();
      int size = m.size();
      buffer.pack<int>(size);
      for (int j = 0; j < size; ++j) {
         buffer.pack<int>(m[j]);
      }

      buffer.incrementSendSize();
   }

   /*
   * Unpack and receive ownership of an Atom.
   */
   void Atom::unpackAtom(Buffer& buffer)
   {
      int i;
      buffer.unpack<int>(i);
      setId(i);
      buffer.unpack<int>(i);
      setTypeId(i);
      buffer.unpack<Vector>(position());
      buffer.unpack<Vector>(velocity());
      unsigned int ui;
      buffer.unpack<unsigned int>(ui);
      plan().setFlags(ui);
      buffer.unpack<unsigned int>(ui);
      groups() = ui;
      if (hasAtomContext_) {
         buffer.unpack<AtomContext>(context());
      }

      // Unpack Mask
      Mask& m = mask();
      m.clear();
      int size;
      buffer.unpack<int>(size);
      for (int j = 0; j < size; ++j) {
         buffer.unpack<int>(i);
         m.append(i);
      }
      assert(m.size() == size);

      buffer.decrementRecvSize();
   }

   /*
   * Return maximum size of packed Atom, in bytes.
   */
   int Atom::packedAtomSize()
   {  
      int size = 0;
      size += 2*sizeof(int);               // id + typeId
      size += 2*sizeof(Vector);            // position + velocity
      size += 2*sizeof(unsigned int);      // plan + groups
      if (hasAtomContext_) {
         size += sizeof(AtomContext);      // context
      }
      size += sizeof(int);                 // mask size
      size += Mask::Capacity*sizeof(int);  // mask ids
      return size;
   }

   // Ghost Exchange

   /*
   * Pack data required for a ghost Atom for sending.
   */
   void Atom::packGhost(Buffer& buffer)
   {
      Vector pos;
      buffer.pack<int>(id());
      buffer.pack<int>(typeId());
      buffer.pack<Vector>(position());
      buffer.pack<unsigned int>(plan().flags());
      if (hasAtomContext_) {
         buffer.pack<AtomContext>(context());
      }
      buffer.incrementSendSize();
   }

   /*
   * Unpack data required for a ghost Atom.
   */
   void Atom::unpackGhost(Buffer& buffer)
   {
      int i;
      buffer.unpack<int>(i);
      setId(i);
      buffer.unpack<int>(i);
      setTypeId(i);
      buffer.unpack<Vector>(position());
      unsigned int ui;
      buffer.unpack<unsigned int>(ui);
      plan().setFlags(ui);
      if (hasAtomContext_) {
         buffer.unpack<AtomContext>(context());
      }
      buffer.decrementRecvSize();
   }

   /*
   * Return size of one packed Ghost  in bytes (static method).
   */
   int Atom::packedGhostSize()
   {  
      int size = 0;
      size += 2*sizeof(int); 
      size += sizeof(Vector); 
      size += sizeof(unsigned int);
      if (hasAtomContext_) {
         size += sizeof(AtomContext);
      }
      return size;
   }

   // Ghost Position Updates

   /*
   * Pack updates ghost position.
   */
   void Atom::packUpdate(Buffer& buffer)
   {
      buffer.pack<Vector>(position());
      buffer.incrementSendSize();
   }

   /*
   * Pack updated ghost position.
   */
   void Atom::unpackUpdate(Buffer& buffer)
   {
      buffer.unpack<Vector>(position());
      buffer.decrementRecvSize();
   }

   // Force Updates

   /*
   * Pack ghost force.
   */
   void Atom::packForce(Buffer& buffer)
   {
      buffer.pack<Vector>(force());
      buffer.incrementSendSize();
   }

   /*
   * Unpack data ghost Atom force, and add to on this processor.
   */
   void Atom::unpackForce(Buffer& buffer)
   {
      Vector f;
      buffer.unpack<Vector>(f);
      force() += f;
      buffer.decrementRecvSize();
   }
   #endif // ifdef UTIL_MPI

   /*
   * Copy ghost atom data from sendAtom to this Atom.
   */
   void Atom::copyLocalGhost(const Atom& sendAtom)
   {
      setId(sendAtom.id());
      setTypeId(sendAtom.typeId());
      plan().setFlags(sendAtom.plan().flags());
      position() = sendAtom.position();
      if (hasAtomContext_) {
         context() = sendAtom.context();
      }
   }

   /*
   * Copy update position of local ghost from sendAtom to this.
   */
   void Atom::copyLocalUpdate(const Atom& sendAtom)
   {  position_ = sendAtom.position(); }

}
