/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "NphIntegrator.h"
#include <mcMd/mdSimulation/MdSystem.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/potentials/pair/MdPairPotential.h>
#include <util/boundary/Boundary.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/ensembles/EnergyEnsemble.h>
#include <util/ensembles/BoundaryEnsemble.h>
#include <util/space/Vector.h>
#include <util/archives/Serializable_includes.h>
#include <util/crystal/LatticeSystem.h>

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor.   
   */
   NphIntegrator::NphIntegrator(MdSystem& system)
    : MdIntegrator(system),
      W_(0.0)
   {
      setClassName("NphIntegrator"); 

      // Start with zero barostat momentum
      eta_.zero();
   }

   /* 
   * Destructor.   
   */
   NphIntegrator::~NphIntegrator() 
   {}

   /* 
   * Read parameter and configuration files, initialize system.
   */
   void NphIntegrator::readParameters(std::istream &in) 
   {
      read<double>(in, "dt", dt_);
      read<double>(in, "W", W_);
      read<LatticeSystem>(in, "mode", mode_);

      int nAtomType = simulation().nAtomType();
      prefactors_.allocate(nAtomType);

   }

   /**
   * Load the internal state to an archive.
   */
   void NphIntegrator::loadParameters(Serializable::IArchive& ar)
   {  
      loadParameter<double>(ar, "dt", dt_);
      loadParameter<double>(ar, "W", W_);
      loadParameter<LatticeSystem>(ar, "mode", mode_);
      ar & eta_;
      ar & currP_;
      ar & prefactors_;
   }

   /*
   * Save the internal state to an archive.
   */
   void NphIntegrator::save(Serializable::OArchive& ar)
   {  
      ar & dt_;
      ar & W_;
      serializeEnum(ar, mode_);
      ar & eta_;
      ar & currP_;
      ar & prefactors_;
   }

   /* 
   * Initialize constants.
   */
   void NphIntegrator::setup() 
   {
      double mass;
      int nAtomType = prefactors_.capacity();
      for (int i = 0; i < nAtomType; ++i) {
         mass = simulation().atomType(i).mass();
         prefactors_[i] = 0.5*dt_/mass;
      }
      system().positionSignal().notify();
      system().velocitySignal().notify();
   }

   void NphIntegrator::step() 
   {
      System::MoleculeIterator molIter;
      double prefactor;
      int    iSpecies, nSpecies;

      nSpecies = simulation().nSpecies();

      /* Perform the first half step of the explicitly reversible NPH integration scheme.
  
         This follows from operator factorization
  
         - orthorhombic case:
  
           1) eta_alpha(t+dt/2) = eta_alpha(t) + dt/2 * V/L_alpha * ( P_{alpha,alpha}(t) - P_ext)
           2) v' = v(t) + (1/2)a(t)*dt
           3) L(t+dt/2) = L(t) + dt*eta(t+dt/2)/(2*W)
           4) r'_alpha = r(t) + v'*dt* (L_{alpha}^2(t)/L_{alpha}^2(t+dt/2))
           5a) L(t+dt) = L(t+dt/2) + dt*eta(t+dt/2)/2/W
           5b) r_alpha(t+dt) = L_alpha(t+dt)/L_alpha(t)*r'_alpha
           5c) v''_alpha = L_alpha(t)/L_alpha(t+dt) * v'_alpha
  
         alpha denotes a cartesian index.
  
         - isotropic case:
  
           only eta_x := eta is used and instead of step 1) we have
  
           1') eta(t+dt/2) = eta(t) + dt/2 * (1/3*Tr P(t) - P_ext)
  
           furthermore, in step 3) and 5a) L is replaced with V=L^3
  
         - tetragonal case:
  
           Lx := L_perp, Ly = Lz := L_par
           eta_x := eta_perp, etay := eta_par, etaz unused
  
           instead of step 1) we have
  
           1'a) eta_perp(t+dt/2) = eta_perp + dt/2 * V/L_perp * ( P_xx(t) - P_ext)
           1'b) eta_par(t+dt/2) = eta_par + dt/2 * V/L_par * ( P_yy(t) + P_zz(t) - 2*P_ext)
  
           steps 3) and 5a) are split into two sub-steps
  
           L_perp(i+1) = L_perp(i) + dt/(2*W)*eta_perp
           L_par(i+1) = L_par(i) + dt/(4*W)*eta_par
      */

      // obtain current stress tensor (diagonal components)
      system().computeStress(currP_);

      // obtain current box lengths
      Vector lengths  = system().boundary().lengths();
      double volume = lengths[0]*lengths[1]*lengths[2];

      double extP = system().boundaryEnsemble().pressure();

      // Advance eta(t)->eta(t+dt/2) (step one)
      if (mode_ == Orthorhombic) {
         double Vdthalf = 1.0/2.0 * volume * dt_;
         eta_[0] += Vdthalf/lengths[0] * (currP_[0] - extP);
         eta_[1] += Vdthalf/lengths[1] * (currP_[1] - extP);
         eta_[2] += Vdthalf/lengths[2] * (currP_[2] - extP);
      } else if (mode_ == Tetragonal) { 
         double Vdthalf = 1.0/2.0* volume * dt_; 
         eta_[0] += Vdthalf/lengths[0]*(currP_[0] - extP);
         eta_[1] += Vdthalf/lengths[1]*(currP_[1] + currP_[2] - 2.0*extP);
      } else if (mode_ == Cubic) {
         eta_[0] += 1.0/2.0*dt_*(1.0/3.0*(currP_[0]+currP_[1]+currP_[2]) - extP);
      }

      // Update the box length L(t) -> L(t+dt/2) (step three)
      // (Since we still keep accelerations a(t) computed for length L_alpha(t) in memory,
      // needed in step two, we can exchange the order of the two steps.)
      // Also pre-calculate L(t+dt) (step 5a, only depends on eta(t) of step one)

      Vector lengthsOld = lengths;
      Vector lengthsFinal;
      double volumeFinal = 0.0;

      double dthalfoverW = 1.0/2.0*dt_/W_;

      if (mode_ == Orthorhombic) {
         lengths[0] += dthalfoverW*eta_[0];
         lengths[1] += dthalfoverW*eta_[1];
         lengths[2] += dthalfoverW*eta_[2];
         lengthsFinal[0] = lengths[0] + dthalfoverW*eta_[0];
         lengthsFinal[1] = lengths[1] + dthalfoverW*eta_[1];
         lengthsFinal[2] = lengths[2] + dthalfoverW*eta_[2];
         volumeFinal = lengthsFinal[0]*lengthsFinal[1]*lengthsFinal[2];
      } else if (mode_ == Tetragonal) {
         lengths[0] += dthalfoverW*eta_[0];
         lengths[1] += (1.0/2.0)*dthalfoverW*eta_[1];
         lengths[2] = lengths[1];
         lengthsFinal[0] = lengths[0] + dthalfoverW*eta_[0];
         lengthsFinal[1] = lengths[1] + (1.0/2.0)*dthalfoverW*eta_[1];
         lengthsFinal[2] = lengthsFinal[1];
         volumeFinal = lengthsFinal[0]*lengthsFinal[1]*lengthsFinal[2];
      } else if (mode_ == Cubic) {
         volume += dthalfoverW*eta_[0];
         lengths[0] = pow(volume,1.0/3.0); // Lx = Ly = Lz = V^(1/3)
         lengths[1] = lengths[0];
         lengths[2] = lengths[0];
         volumeFinal = volume + dthalfoverW*eta_[0];
         lengthsFinal[0] = pow(volumeFinal,1./3.);
         lengthsFinal[1] = lengthsFinal[0];
         lengthsFinal[2] = lengthsFinal[0];
      }
         
      // Update simulation box 
      system().boundary().setOrthorhombic(lengthsFinal); 
      
      // update particles
      Atom* atomPtr;
      int   ia;
      Vector vtmp, dr, rtmp, dv;
      Vector lengthsFac, dtLengthsFac2;

      for (int i=0; i<3; i++) {
         lengthsFac[i] = lengthsFinal[i]/lengthsOld[i];
         dtLengthsFac2[i] = dt_ * lengthsOld[i]*lengthsOld[i]/lengths[i]/lengths[i];
      }

      for (iSpecies=0; iSpecies < nSpecies; ++iSpecies) {
         system().begin(iSpecies, molIter); 
         for ( ; molIter.notEnd(); ++molIter) {
            for (ia=0; ia < molIter->nAtom(); ++ia) {
               atomPtr = &molIter->atom(ia);
               prefactor = prefactors_[atomPtr->typeId()];

               dv.multiply(atomPtr->force(), prefactor);
               vtmp.add(atomPtr->velocity(),dv);

               dr[0] = vtmp[0] * dtLengthsFac2[0];
               dr[1] = vtmp[1] * dtLengthsFac2[1];
               dr[2] = vtmp[2] * dtLengthsFac2[2];

               rtmp.add(atomPtr->position(), dr);

               rtmp[0] *= lengthsFac[0];
               rtmp[1] *= lengthsFac[1];
               rtmp[2] *= lengthsFac[2];

               atomPtr->position() = rtmp;

               vtmp[0] /= lengthsFac[0];
               vtmp[1] /= lengthsFac[1];
               vtmp[2] /= lengthsFac[2];

               atomPtr->velocity() = vtmp;
            }
         }
      }
      system().positionSignal().notify();
      system().velocitySignal().notify();

      #ifndef SIMP_NOPAIR
      if (!system().pairPotential().isPairListCurrent()) {
         system().pairPotential().buildPairList();
      } 
      #endif

      system().calculateForces();

      /* Second step of the explicitly reversible integrator consists of the following sub-steps
      *
      *   6) v(t+dt) = v'' + 1/2 * a(t+dt)*dt
      *   7) eta(t+dt/2) -> eta(t+dt)
      */

      // Update velocities
      // v(t+dt) = v'' + 1/2 * a(t+dt)*dt
      for (iSpecies=0; iSpecies < nSpecies; ++iSpecies) {
         system().begin(iSpecies, molIter); 
         for ( ; molIter.notEnd(); ++molIter)
         {
            for (ia=0; ia < molIter->nAtom(); ++ia) {
               atomPtr = &molIter->atom(ia);
               prefactor = prefactors_[atomPtr->typeId()];
               dv.multiply(atomPtr->force(), prefactor);
               atomPtr->velocity() += dv;
            }
         }
      }
      system().velocitySignal().notify();

      // now compute pressure tensor with updated virial and velocities
      system().computeStress(currP_);

      //  Advance eta(t+dt/2) -> eta(t+dt)
      if (mode_ == Orthorhombic) {
         double Vdthalf = 1.0/2.0 * volumeFinal * dt_;
         eta_[0] += Vdthalf/lengthsFinal[0] * (currP_[0] - extP);
         eta_[1] += Vdthalf/lengthsFinal[1] * (currP_[1] - extP);
         eta_[2] += Vdthalf/lengthsFinal[2] * (currP_[2] - extP);
      } else if (mode_ == Tetragonal) {
         double Vdthalf = 1.0/2.0 * volumeFinal * dt_;
         eta_[0] += Vdthalf/lengthsFinal[0]*(currP_[0] - extP);
         eta_[1] += Vdthalf/lengthsFinal[1]*(currP_[1] + currP_[2] - 2.0*extP);
      } else if (mode_ == Cubic) {
         eta_[0] += 1.0/2.0*dt_*(1.0/3.0*(currP_[0]+currP_[1]+currP_[2]) - extP);
      }
 
      #ifndef SIMP_NOPAIR
      if (!system().pairPotential().isPairListCurrent()) {
         system().pairPotential().buildPairList();
      }
      #endif
     
      /*
      // debug output
      double P;
      system().computeStress(P);
      std::cout << system().boundary() << std::endl;
      std::cout << "P = " << P << " ";
      std::cout << "Ekin = " << system().kineticEnergy() << " ";
      std::cout << "Epot = " << system().potentialEnergy() << " ";
      std::cout << "PV = " << extP * system().boundary().volume() << " ";
      std::cout << "Ebaro = " << barostatEnergy() << " "; 
      std::cout << "H = " << system().potentialEnergy() + system().kineticEnergy() +
                    extP * system().boundary().volume() + barostatEnergy() << std::endl;
      */
   }

   /* 
   * Get barostat energy.
   */
   double NphIntegrator::barostatEnergy()
   {
      double barostat_energy = 0.0;
      if (mode_ == Orthorhombic)
         barostat_energy = (eta_[0]*eta_[0]+eta_[1]*eta_[1]+eta_[2]*eta_[2]) / 2.0 / W_;
      else if (mode_ == Tetragonal)
         barostat_energy = (eta_[0]*eta_[0]+eta_[1]*eta_[1]/2.0) / 2.0 / W_;
      else if (mode_ == Cubic)
         barostat_energy = eta_[0]*eta_[0]/2.0/W_;

      return barostat_energy;
   }

   /* 
   * Get barostat mass.
   */
   double NphIntegrator::barostatMass() const
   {  return W_; }
 
   /* 
   * Get mode of integrator.
   */
   LatticeSystem NphIntegrator::mode() const
   {  return mode_; }

   /* 
   * Get momentum of barostat.
   */
   //Vector& NphIntegrator::eta()
   //{  return eta_; }

   /* 
   * Set momentum component of barostat.
   */
   void NphIntegrator::setEta(unsigned int index, double eta)
   {  
      eta_[index] = eta; 
   }
   
}
