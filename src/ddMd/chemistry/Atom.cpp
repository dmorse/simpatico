#ifndef DDMD_ATOM_CPP
#define DDMD_ATOM_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
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
   * Constructor (private, used by AtomArray).
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
      setIsGhost(other.isGhost());
      force_ = other.force_;
      velocity() = other.velocity();
      #ifdef DDMD_MOLECULES
      context() = other.context();
      #endif
      mask() = other.mask();
      plan() = other.plan();
      setId(other.id());
      return *this;
   }

   /*
   * Reset integer members to null values.
   */
   void Atom::clear()
   {
      typeId_ = -1;
      setIsGhost(false);
      mask().clear();
      plan().clearFlags();
      setId(-1);
   }

   #ifdef UTIL_MPI
   /*
   * Pack a local Atom for exchange of ownership.
   */
   void Atom::packAtom(Buffer& buffer)
   {
      buffer.pack<int>(id());
      buffer.pack<int>(typeId());
      buffer.pack<Vector>(position());
      buffer.pack<Vector>(velocity());
      #ifdef DDMD_MOLECULES
      buffer.pack<AtomContext>(context());
      #endif
      buffer.pack<unsigned int>(plan().flags());

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
   * Receive ownership of an Atom.
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
      #ifdef DDMD_MOLECULES
      buffer.unpack<AtomContext>(context());
      #endif
      unsigned int ui;
      buffer.unpack<unsigned int>(ui);
      plan().setFlags(ui);

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
   * Return size of packed Atom, in bytes.
   */
   int Atom::packedAtomSize()
   {  
      int size = 0;
      size += 2*sizeof(int);
      #ifdef DDMD_MOLECULES 
      size += sizeof(AtomContext);
      #endif 
      size += 2*sizeof(Vector); 
      size += sizeof(unsigned int);
      size += sizeof(int);
      size += Mask::Capacity*sizeof(int); 
      return size;
   }

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
      return size;
   }

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
   #endif

}
#endif
