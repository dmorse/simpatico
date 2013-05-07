#ifndef DDMD_INTEGRATOR_H
#define DDMD_INTEGRATOR_H

#include <util/param/ParamComposite.h>          // base class
#include <ddMd/simulation/SimulationAccess.h>   // base class
#include <ddMd/misc/DdTimer.h>                  // member

#include <iostream>

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace DdMd
{

   using namespace Util;

   class Simulation;

   /**
   * An Integrator numerically integrates the equations of motion.
   *
   * \ingroup DdMd_Integrator_Module
   */
   class Integrator : public ParamComposite, public SimulationAccess
   {

   public:

      /**
      * Constructor.
      */
      Integrator(Simulation& simulation);

      /**
      * Destructor.
      */
      ~Integrator();

      /**
      * Read saveInterval and saveFileName.
      *
      * \param in input parameter stream
      */   
      void readParameters(std::istream& in);

      /**
      * Load saveInterval and saveFileName from restart archive.
      *
      * \param ar input archive
      */   
      void loadParameters(Serializable::IArchive& ar);

      /**
      * Save saveInterval and saveFileName from restart archive.
      *
      * \param ar input archive
      */   
      void save(Serializable::OArchive& ar);

      /**
      * Run a simulation of iStep steps.
      *
      * \param nStep number of steps to run.
      */
      virtual void run(int nStep) = 0;

      /**
      * Set integrator to initial state and clears all statistics.
      *
      * This method resset iStep = 0, calls initDynamicalState(), calls the
      * DiagnosticManger::clear(), clears the internal Timer and all timing 
      * statistics, clear the additional timing statistics maintained by
      * the Exchanger class, and clears the memory usage statistics for the 
      * Buffer, PairList, and storage classes. 
      */
      virtual void clear();

      /*
      * Reduce timing statistics data from all domain processors.
      *
      * This method must be called simultaneously on all processors. It should
      * always be called immediately before calling outputStatistics on the
      * master processor. 
      *
      * On return, the value of each of the time intervals stored by the
      * internal DdMd::Timer is replaced by its average over all processors.
      * It does not reset the timer to zero, and may thus be called more 
      * than once while contining to accumulate statistics.
      */
      void computeStatistics();

      /**
      * Output timing statistics from a run.
      *
      * This method may only be called on the domain master processor.
      *
      * \param out output stream to which to write timing statistics.
      */
      virtual void outputStatistics(std::ostream& out);

      /**
      * Get average time per processor of previous run.
      */
      double time() const;

      /**
      * Get current time step index.
      */
      int iStep() const;

   protected:

      /// Timestamps for loop timing.
      enum TimeId {DIAGNOSTIC, INTEGRATE1, CHECK, ALLREDUCE, TRANSFORM_F, 
                   EXCHANGE, CELLLIST, TRANSFORM_R, PAIRLIST, UPDATE, 
                   PAIR_FORCE, BOND_FORCE, ANGLE_FORCE, DIHEDRAL_FORCE,
                   INTEGRATE2, NTime};

      /**
      * Set any internal dynamical to default initial values.
      *
      * This method should be called before the main loop the first time the 
      * integrator is used, within the setup() method, and should be called 
      * by the clear() method.
      */
      virtual void initDynamicalState();

      /**
      * Setup integrator just before integration.
      *
      * This method is always called within the run method before the main loop.
      * It should not set arbitrary default values for any independent internal 
      * state variables that cannot be calculated from the system configuration.
      */
      virtual void setup() = 0;

      /**
      * Mark the integrator as having been setup at least once.
      *
      * Must be called from within the setup method.
      */
      void setIsSetup();

      /**
      * Has the setup() method been called at least once previously?
      */
      bool isSetup() const;

      /**
      * Setup state of atoms just before integration.
      *
      * Exchange atoms, build PairList and compute Forces.
      * Should be called in all subclass setup methods.
      */
      void setupAtoms();

      /**
      * Compute forces for all local atoms, with timing.
      *
      * Identical to Simulation::forceCompute(), with added timing.
      * Upon return, forces are correct for all local atoms. Values
      * of the forces on ghost atoms are undefined. Executes reverse
      * communication if needed, and emits Simulation::forceSignal().
      */
      void computeForces();

      /**
      * Compute forces for all local atoms and virial, with timing.
      *
      * Identical to Simulation::computeForcesAndVirial(), with timing.
      * Upon return, forces are correct for all local atoms. Values
      * of the forces on ghost atoms are undefined. Executes reverse
      * communication if needed, and emits Simulation::forceSignal().
      */
      void computeForcesAndVirial();

      /**
      * Determine whether an atom exchange and reneighboring is needed.
      *
      * \param skin Verlet list skin length
      * \return true iff exchange is needed
      */
      bool isExchangeNeeded(double skin);

      /**
      * Get restart file base name. 
      */
      const std::string& saveFileName() const;

      /**
      * Get interval for writing a restart file.
      */
      int saveInterval() const;

      /*
      * Return the timer by reference.
      */
      DdTimer& timer();

      /// Current step number.
      int   iStep_;

   private:

      // Performance timer
      DdTimer timer_;

      // Has setup been called at least once?
      bool isSetup_;

      /// Restart output file name
      std::string saveFileName_;

      /// Interval for writing restart files (no output if 0)
      int saveInterval_;

   };

   /*
   * Get current step index.
   */
   inline int Integrator::iStep() const
   {  return iStep_ ; }

   /*
   * Has setup been called at least once?
   */
   inline bool Integrator::isSetup() const
   {  return isSetup_ ; }

   /*
   * Mark as having been setup at least once.
   */
   inline void Integrator::setIsSetup()
   {  isSetup_ = true; }

   /*
   * Set any internal dynamical state variables to default initial values.
   *  
   * Empty default implementation. This virtual method is called by clear().
   */
   inline void Integrator::initDynamicalState(){}

   /*
   * Return the timer by reference.
   */
   inline DdTimer& Integrator::timer()
   {  return timer_; }

   /*
   * Get restart file base name. 
   */
   inline const std::string& Integrator::saveFileName() const
   { return saveFileName_; }

   /*
   * Get interval for writing a restart file.
   */
   inline int Integrator::saveInterval() const
   { return saveInterval_; }

}
#endif
