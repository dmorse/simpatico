/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "NphIntegrator.h"
#include <ddMd/simulation/Simulation.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/communicate/Exchanger.h>
#include <ddMd/potentials/pair/PairPotential.h>
#include <util/space/Vector.h>
#include <util/ensembles/EnergyEnsemble.h>
#include <util/ensembles/BoundaryEnsemble.h>
#include <util/mpi/MpiLoader.h>
#include <util/format/Dbl.h>
#include <util/global.h>

#include <iostream>

namespace DdMd
{
   using namespace Util;

   /*
   * Constructor.
   */
   NphIntegrator::NphIntegrator(Simulation& simulation)
    : TwoStepIntegrator(simulation)
   { setClassName("NphIntegrator"); }

   /*
   * Destructor.
   */
   NphIntegrator::~NphIntegrator()
   {}

   /*
   * Read parameters and initialize.
   */
   void NphIntegrator::readParameters(std::istream& in)
   {
      read<double>(in, "dt", dt_);
      read<double>(in, "W", W_);
      read<LatticeSystem>(in, "mode", mode_);
      Integrator::readParameters(in);

      // Allocate memory
      int nAtomType = simulation().nAtomType();
      if (!prefactors_.isAllocated()) {
         prefactors_.allocate(nAtomType);
      }
   }

   /*
   * Load state from restart archive.
   */
   void NphIntegrator::loadParameters(Serializable::IArchive& ar)
   {
      loadParameter<double>(ar, "dt", dt_);
      loadParameter<double>(ar, "W", W_);
      loadParameter<LatticeSystem>(ar, "mode", mode_);
      Integrator::loadParameters(ar);

      MpiLoader<Serializable::IArchive> loader(*this, ar);
      loader.load(nu_);

      // Allocate memory
      int nAtomType = simulation().nAtomType();
      if (!prefactors_.isAllocated()) {
         prefactors_.allocate(nAtomType);
      }
   }

   /*
   * Save state to output archive.
   */
   void NphIntegrator::save(Serializable::OArchive& ar)
   {
      ar << dt_;
      ar << W_;
      ar << mode_;
      Integrator::save(ar);
      ar << nu_;
   }

   /*
   * Initialize extended ensemble state variables.
   */
   void NphIntegrator::initDynamicalState()
   {  nu_ = Vector(0.0,0.0,0.0); }

   void NphIntegrator::setup()
   {
      // Initialize state and clear statistics on first usage.
      if (!isSetup()) {
         clear();
         setIsSetup();
      }

      // Exchange atoms, build pair list, compute forces.
      setupAtoms();

      // Set prefactors for acceleration
      double dtHalf = 0.5*dt_;
      double mass;
      int nAtomType = prefactors_.capacity();
      for (int i = 0; i < nAtomType; ++i) {
         mass = simulation().atomType(i).mass();
         prefactors_[i] = dtHalf/mass;
      }

      simulation().computeKineticEnergy();
      #ifdef UTIL_MPI
      atomStorage().computeNAtomTotal(domain().communicator());
      #endif
      if (domain().isMaster()) {
         P_target_ = simulation().boundaryEnsemble().pressure();
         ndof_ = atomStorage().nAtomTotal()*3;
      }
      #ifdef UTIL_MPI
      bcast(domain().communicator(), ndof_,0);
      #endif
   }

   /*
   *  First half of two-step velocity-Verlet integrator. 
   */
   void NphIntegrator::integrateStep1()
   {
      Vector dv;
      AtomIterator atomIter;

      Simulation& sys = simulation();
      sys.computeVirialStress();
      sys.computeKineticStress();
      sys.computeKineticEnergy();

      if (sys.domain().isMaster()) {
         P_target_ = simulation().boundaryEnsemble().pressure();
         T_kinetic_ = sys.kineticEnergy()*2.0/ndof_;
         Tensor stress = sys.virialStress();
         stress += sys.kineticStress();

         P_curr_diag_ = Vector(stress(0, 0), stress(1,1), stress(2,2));
         double P_curr = (1.0/3.0)*stress.trace();

         double mtk_term = (1.0/2.0)*dt_*T_kinetic_/W_;

         // Advance barostat
         double V = sys.boundary().volume();
         if (mode_ == Cubic) {
            nu_[0] += (1.0/2.0)*dt_*V/W_*(P_curr - P_target_) + mtk_term;
            nu_[1] = nu_[2] = nu_[0];
         } else if (mode_ == Tetragonal) {
            nu_[0] += (1.0/2.0)*dt_*V/W_*(P_curr_diag_[0] - P_target_) + mtk_term;
            nu_[1] += (1.0/2.0)*dt_*V/W_*((1.0/2.0)*(P_curr_diag_[1]+P_curr_diag_[2]) - P_target_) + mtk_term;
            nu_[2] = nu_[1];
         } else if (mode_  == Orthorhombic) {
            nu_[0] += (1.0/2.0)*dt_*V/W_*(P_curr_diag_[0] - P_target_) + mtk_term;
            nu_[1] += (1.0/2.0)*dt_*V/W_*(P_curr_diag_[1] - P_target_) + mtk_term;
            nu_[2] += (1.0/2.0)*dt_*V/W_*(P_curr_diag_[2] - P_target_) + mtk_term;
         }

      }

      #ifdef UTIL_MPI
      bcast(domain().communicator(), nu_,0);
      #endif

      // Precompute loop invariant quantities
      mtk_term_2_ = (nu_[0]+nu_[1]+nu_[2])/ndof_;
      Vector v_fac = Vector((1.0/4.0)*(nu_[0]+mtk_term_2_),
                            (1.0/4.0)*(nu_[1]+mtk_term_2_),
                            (1.0/4.0)*(nu_[2]+mtk_term_2_));
      exp_v_fac_ = Vector(exp(-v_fac[0]*dt_),
                          exp(-v_fac[1]*dt_),
                          exp(-v_fac[2]*dt_));
      Vector exp_v_fac_2 = Vector(exp(-(2.0*v_fac[0])*dt_),
                                  exp(-(2.0*v_fac[1])*dt_),
                                  exp(-(2.0*v_fac[2])*dt_));
      Vector r_fac = Vector((1.0/2.0)*nu_[0],
                            (1.0/2.0)*nu_[1],
                            (1.0/2.0)*nu_[2]);
      Vector exp_r_fac = Vector(exp(r_fac[0]*dt_),
                                exp(r_fac[1]*dt_),
                                exp(r_fac[2]*dt_));

      // Coefficients of sinh(x)/x = a_0 + a_2*x^2 + a_4*x^4 + a_6*x^6 + a_8*x^8 + a_10*x^10
      const double a[] = {double(1.0), double(1.0/6.0), double(1.0/120.0), 
                          double(1.0/5040.0), double(1.0/362880.0), double(1.0/39916800.0)};

      Vector arg_v = Vector(v_fac[0]*dt_, v_fac[1]*dt_, v_fac[2]*dt_);
      Vector arg_r = Vector(r_fac[0]*dt_, r_fac[1]*dt_, r_fac[2]*dt_);

      sinhx_fac_v_ = Vector(0.0,0.0,0.0);
      Vector sinhx_fac_r = Vector(0.0,0.0,0.0);

      Vector term_v = Vector(1.0,1.0,1.0);
      Vector term_r = Vector(1.0,1.0,1.0);

      for (unsigned int i = 0; i < 6; i++) {
         sinhx_fac_v_ += Vector(a[i]*term_v[0],
                                a[i]*term_v[1],
                                a[i]*term_v[2]);
         sinhx_fac_r += Vector(a[i]*term_r[0],
                               a[i]*term_r[1],
                               a[i]*term_r[2]);
         term_v = Vector(term_v[0] * arg_v[0] * arg_v[0],
                         term_v[1] * arg_v[1] * arg_v[1],
                         term_v[2] * arg_v[2] * arg_v[2]);
         term_r = Vector(term_r[0] * arg_r[0] * arg_r[0],
                         term_r[1] * arg_r[1] * arg_r[1],
                         term_r[2] * arg_r[2] * arg_r[2]);
      }

      // 1st half of NPH
      Vector vtmp;
      double prefactor; // = 0.5*dt/mass
      atomStorage().begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         prefactor = prefactors_[atomIter->typeId()];

         dv.multiply(atomIter->force(), prefactor);
         dv[0] = dv[0] * exp_v_fac_[0]*sinhx_fac_v_[0];
         dv[1] = dv[1] * exp_v_fac_[1]*sinhx_fac_v_[1];
         dv[2] = dv[2] * exp_v_fac_[2]*sinhx_fac_v_[2];

         Vector& v = atomIter->velocity();
         v[0] = v[0] * exp_v_fac_2[0] + dv[0];
         v[1] = v[1] * exp_v_fac_2[1] + dv[1];
         v[2] = v[2] * exp_v_fac_2[2] + dv[2];

         vtmp[0] = v[0]*exp_r_fac[0] *sinhx_fac_r[0];
         vtmp[1] = v[1]*exp_r_fac[1] *sinhx_fac_r[1];
         vtmp[2] = v[2]*exp_r_fac[2] *sinhx_fac_r[2];

         Vector& r = atomIter->position();
         r[0] = r[0]*exp_r_fac[0]*exp_r_fac[0] + vtmp[0]*dt_;
         r[1] = r[1]*exp_r_fac[1]*exp_r_fac[1] + vtmp[1]*dt_;
         r[2] = r[2]*exp_r_fac[2]*exp_r_fac[2] + vtmp[2]*dt_;
      }

      // Advance box lengths
      Vector box_len_scale = Vector(exp(nu_[0]*dt_),
                                    exp(nu_[1]*dt_),
                                    exp(nu_[2]*dt_));

      Vector L = sys.boundary().lengths();
      L[0] *= box_len_scale[0];
      L[1] *= box_len_scale[1];
      L[2] *= box_len_scale[2];

      // Update box lengths
      sys.boundary().setOrthorhombic(L);
   }

   /*
   *  Second half of two-step velocity-Verlet integrator. 
   */
   void NphIntegrator::integrateStep2()
   {
      Vector dv;
      double prefactor; // = 0.5*dt/mass
      AtomIterator atomIter;

      Vector v_fac_2 = Vector((1.0/2.0)*(nu_[0]+mtk_term_2_),
                              (1.0/2.0)*(nu_[1]+mtk_term_2_),
                              (1.0/2.0)*(nu_[2]+mtk_term_2_));
      Vector exp_v_fac_2 = Vector(exp(-v_fac_2[0]*dt_),
                                 exp(-v_fac_2[1]*dt_),
                                 exp(-v_fac_2[2]*dt_));

      // 2nd half of NPH
      atomStorage().begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         prefactor = prefactors_[atomIter->typeId()];
         dv.multiply(atomIter->force(), prefactor);
         dv[0] = dv[0] * exp_v_fac_[0]*sinhx_fac_v_[0];
         dv[1] = dv[1] * exp_v_fac_[1]*sinhx_fac_v_[1];
         dv[2] = dv[2] * exp_v_fac_[2]*sinhx_fac_v_[2];

         Vector& v = atomIter->velocity();
         v[0] = v[0] * exp_v_fac_2[0] + dv[0];
         v[1] = v[1] * exp_v_fac_2[1] + dv[1];
         v[2] = v[2] * exp_v_fac_2[2] + dv[2];
      }

      Simulation& sys = simulation();
      sys.velocitySignal().notify();
      sys.computeKineticStress();  
      sys.computeKineticEnergy();  
      sys.computeVirialStress(); 

      // Advance barostat
      if (sys.domain().isMaster()) {
         T_kinetic_ = sys.kineticEnergy()*2.0/ndof_;
         Tensor stress = sys.virialStress();
         stress += sys.kineticStress();

         P_curr_diag_ = Vector(stress(0,0), stress(1,1), stress(2,2));
         double P_curr = (1.0/3.0)*stress.trace();

         double mtk_term = (1.0/2.0)*dt_*T_kinetic_/W_;

         double V = sys.boundary().volume();
         if (mode_ == Cubic) {
            nu_[0] += (1.0/2.0)*dt_*V/W_*(P_curr - P_target_) + mtk_term;
            nu_[1] = nu_[2] = nu_[0];
         } else if (mode_ == Tetragonal) {
            nu_[0] += (1.0/2.0)*dt_*V/W_*(P_curr_diag_[0] - P_target_) + mtk_term;
            nu_[1] += (1.0/2.0)*dt_*V/W_*((1.0/2.0)*(P_curr_diag_[1]+P_curr_diag_[2]) - P_target_) + mtk_term;
            nu_[2] = nu_[1];
         } else if (mode_  == Orthorhombic) {
            nu_[0] += (1.0/2.0)*dt_*V/W_*(P_curr_diag_[0] - P_target_) + mtk_term;
            nu_[1] += (1.0/2.0)*dt_*V/W_*(P_curr_diag_[1] - P_target_) + mtk_term;
            nu_[2] += (1.0/2.0)*dt_*V/W_*(P_curr_diag_[2] - P_target_) + mtk_term;
         }
      }

      #if 0
      // output conserved quantity
      sys.computePotentialEnergies();
      if (sys.domain().isMaster()) {
         std::ofstream file;
         file.open("NPH.log", std::ios::out | std::ios::app);
         double V = sys.boundary().volume();
         double barostat_energy = W_*(nu_[0]*nu_[0]+ nu_[1]*nu_[1] + nu_[2]*nu_[2]);
         double pe = sys.potentialEnergy();
         double ke = sys.kineticEnergy();
         file << Dbl(V,20)
              << Dbl(pe,20)
              << Dbl(ke,20)
              << Dbl(barostat_energy,20)
              << std::endl;
         file.close();
      }
      #endif

      #ifdef UTIL_MPI
      bcast(domain().communicator(), nu_,0);
      #endif

  }

}
