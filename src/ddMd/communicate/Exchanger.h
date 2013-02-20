#ifndef DDMD_EXCHANGER_H
#define DDMD_EXCHANGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/chemistry/Atom.h>
#include <ddMd/chemistry/Group.h>
#include <ddMd/misc/DdTimer.h>
#include <util/boundary/Boundary.h>
#include <util/containers/FMatrix.h>
#include <util/containers/APArray.h>


namespace DdMd
{

   class Domain;
   class AtomStorage;
   class Buffer;
   template <int N> class GroupStorage;

   using namespace Util;

   /**
   * Class for exchanging Atoms, Ghosts and Groups between processors.
   *
   * \ingroup DdMd_Communicate_Module
   */
   class Exchanger
   {

   public:

      /**
      * Constructor.
      */
      Exchanger();

      /**
      * Destructor.
      */
      ~Exchanger();

      /**
      * Set pointers to associated objects.
      */
      void associate(const Domain& domain, 
                     const Boundary& boundary, 
                     AtomStorage& atomStorage, 
                     GroupStorage<2>& bondStorage, 
                     #ifdef INTER_ANGLE
                     GroupStorage<3>& angleStorage, 
                     #endif
                     #ifdef INTER_DIHEDRAL
                     GroupStorage<4>& dihedralStorage, 
                     #endif
                     Buffer& buffer);

      /**
      * Allocate all required memory.
      *
      * Must be called after Buffer is initialized.
      */
      void allocate();

      /**
      * Set width of slab for ghosts.
      *
      * \param pairCutoff cutoff radius for pair list (potential + skin).
      */
      void setPairCutoff(double pairCutoff);

      /**
      * Exchange local atoms and ghosts.
      * 
      * This method should be called just before rebuilding the neighbor
      * list on each processor, to exchange ownership of local atoms and
      * to exchange ghost atoms. The lists of which atoms were sent and
      * received as ghosts by this method are used in subsequent calls
      * to update. 
      */
      void exchange();

      /**
      * Update ghost atom coordinates.
      * 
      * This method should be called every time step for which there is
      * no exhange of atom ownership. It communicates ghost coordinates
      * for the same ghosts as those sent by the most recent call to
      * the exchangeGhosts() methods.
      */
      void update();

      /**
      * Update ghost atom forces.
      * 
      * This method reverse the communication pattern used to communicate
      * ghost atom positions in update() to reverse communicate forces
      * acting on ghost atoms. It should be called iff reverse force
      * communication on every time step for which update() is called.
      */
      void reverseUpdate();

      /**
      * Output statistics.
      */
      void outputStatistics(std::ostream& out, double time, int nStep);

      /**
      * Return internal timer by reference
      */
      DdTimer& timer()
      {  return timer_; }

      /**
      * Enumeration of time stamp identifiers.
      */
      enum timeId {START, ATOM_PLAN, INIT_GROUP_PLAN, CLEAR_GHOSTS,
                   PACK_ATOMS, PACK_GROUPS, REMOVE_ATOMS, REMOVE_GROUPS,
                   SEND_RECV_ATOMS, UNPACK_ATOMS, UNPACK_GROUPS, 
                   FINISH_GROUP_PLAN, INIT_SEND_ARRAYS, PACK_GHOSTS, 
                   SEND_RECV_GHOSTS, UNPACK_GHOSTS, FIND_GROUP_GHOSTS, 
                   PACK_UPDATE, SEND_RECV_UPDATE, UNPACK_UPDATE, 
                   LOCAL_UPDATE, PACK_FORCE, SEND_RECV_FORCE, 
                   UNPACK_FORCE, LOCAL_FORCE, NTime};

   private:

      /**
      * Arrays of pointers to ghosts to be sent in each direction.
      *
      * Element sendArray_(i, j) is a APArray that stores pointers to the
      * ghost atoms whose positions are sent during the send along cartesian 
      * axis i (i=0/x, 1/y or 2/z) in direction j (j=0/lower or 1/higher).
      * These 6 arrays are filled by the exchangeGhosts() method and used 
      * by subsequent calls of the update() method, until the next
      * call of exchangeGhosts().
      */
      FMatrix< APArray<Atom>, Dimension, 2>  sendArray_;

      /**
      * Arrays of pointers to ghosts received in each direction.
      *
      * Element recvArray_(i, j) is a APArray that stores pointers to the 
      * ghost atoms that are received during the send along cartesian axis
      * i (for i=0/x, 1/y or 2/z) in direction j (for j=0/lower or 1/higher).
      * These 6 arrays are filled by the exchangeGhosts() method and used 
      * by subsequent calls of the update() method, until the next
      * call of exchangeGhosts().
      */
      FMatrix< APArray<Atom>, Dimension, 2>  recvArray_;

      #ifdef UTIL_MPI
      /**
      * Array of pointers to atoms that have been packed and sent.
      * 
      * Used to mark missing atoms for subsequent removal.
      */
      APArray<Atom> sentAtoms_;

      /**
      * Array of pointers to empty bonds on this processor.
      * 
      * Used to mark bonds for later removal.
      */
      APArray< Group<2> > emptyBonds_;

      #ifdef INTER_ANGLE
      /**
      * Array of pointers to empty angles on this processor.
      */
      APArray< Group<3> > emptyAngles_;
      #endif

      #ifdef INTER_DIHEDRAL
      /**
      * Array of pointers to empty dihedrals on this processor.
      */
      APArray< Group<4> > emptyDihedrals_;
      #endif

      #endif // UTIL_MPI

      /// Processor boundaries (minima j=0, maxima j=1)
      FMatrix< double, Dimension, 2>  bound_;

      /// Inner boundaries of nonbonded slabs
      FMatrix< double, Dimension, 2>  inner_;

      /// Outer boundaries of nonbonded slabs
      FMatrix< double, Dimension, 2>  outer_;

      /// Elements are 1 if grid dimension > 1, 0 otherwise.
      IntVector multiProcessorDirection_;

      /// Pointer to associated const Boundary object.
      const Boundary*  boundaryPtr_;

      /// Pointer to associated const Domain object.
      const Domain*  domainPtr_;

      /// Pointer to associated AtomStorage object.
      AtomStorage* atomStoragePtr_;

      /// Pointer to associated bond storage object.
      GroupStorage<2>* bondStoragePtr_;

      #ifdef INTER_ANGLE 
      /// Pointer to associated angle storage object.
      GroupStorage<3>* angleStoragePtr_;
      #endif

      #ifdef INTER_DIHEDRAL
      /// Pointer to associated dihedral storage object.
      GroupStorage<4>* dihedralStoragePtr_;
      #endif

      /// Pointer to associated buffer object.
      Buffer*  bufferPtr_;

      /// Cutoff for pair list (potential cutoff + skin).
      double pairCutoff_;

      /// Timer
      DdTimer timer_;

      /**
      * Exchange ownership of local atoms.
      *
      * This method exchanges ownershp of local atoms, and calculates
      * communication plan for ghost atoms, but but does not actually
      * exchange ghost atoms.  
      */
      void exchangeAtoms();

      /**
      * Exchange ghosts.
      *
      * This method exchanges ghosts, and stores lists of which atoms 
      * are sent and received in each direction for use in subsequent
      * calls to update(). It must be called immediately after
      * exchangeAtoms().
      */
      void exchangeGhosts();

      template <int N>
      void initGroupGhostPlan(GroupStorage<N>& storage);

      #ifdef UTIL_MPI
      template <int N>
      void packGroups(int i, int j, GroupStorage<N>& storage, 
                                    APArray< Group<N> >& emptyGroups);

      template <int N>
      void removeEmptyGroups(GroupStorage<N>& storage,
                             APArray< Group<N> >& emptyGroups);

      template <int N>
      void unpackGroups(GroupStorage<N>& storage);
      #endif

      template <int N>
      void finishGroupGhostPlan(GroupStorage<N>& storage);

      template <int N>
      void findGroupGhosts(GroupStorage<N>& storage);

      void stamp(unsigned int timeId) 
      {  timer_.stamp(timeId); }

   };

}
#endif
