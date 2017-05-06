#ifndef MD_EWALD_POTENTIAL_CPP
#define MD_EWALD_POTENTIAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "MdEwaldPotential.h" 
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <util/space/Vector.h>
#include <util/space/Tensor.h>
#include <util/math/Constants.h>
#include <stdlib.h>

#include <util/boundary/Boundary.h>
#include <cmath>
#include <util/containers/Array.h>
#include <mcMd/chemistry/AtomType.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   MdEwaldPotential::MdEwaldPotential(System& system)
    : MdCoulombPotential(),
      ewaldInteraction_(),
      simulationPtr_(&system.simulation()),
      systemPtr_(&system),
      boundaryPtr_(&system.boundary()),
      atomTypesPtr_(&system.simulation().atomTypes())
   {
      //initialize unit tensor.
      const double unitMatrix_[3][3] = { {1,0,0}, {0,1,0}, {0,0,1}};
      Tensor unitTensor_(unitMatrix_);
      // Note: Don't setClassName - using "CoulombPotential" base class name
   }

   /*
   * Destructor (does nothing)
   */
   MdEwaldPotential::~MdEwaldPotential()
   {}

   /*
   * Read parameters and initialize.
   */
   void MdEwaldPotential::readParameters(std::istream& in)
   {
      // Read EwaldInteraction block containing parameters
      bool nextIndent = false;
      addParamComposite(ewaldInteraction_, nextIndent);
      ewaldInteraction_.readParameters(in);

      read<double>(in, "kSpaceCutoff", kSpaceCutoff_);
      kSpaceCutoffSq_ = kSpaceCutoff_*kSpaceCutoff_;

      //Calculate prefactors frequently used in this class.
      double pi = Constants::Pi;
      selfPrefactor_ = ewaldInteraction_.alpha()/(4*sqrt(pi)*pi*ewaldInteraction_.epsilon());
   }

   /*
   * Load internal state from an archive.
   */
   void MdEwaldPotential::loadParameters(Serializable::IArchive &ar)
   {
      bool nextIndent = false;
      addParamComposite(ewaldInteraction_, nextIndent);
      ewaldInteraction_.loadParameters(ar);

      loadParameter<double>(ar, "kSpaceCutoff", kSpaceCutoff_);
      kSpaceCutoffSq_ = kSpaceCutoff_ * kSpaceCutoff_;

      // Calculate selfPrefactor_
      double pi = Constants::Pi;
      selfPrefactor_ = ewaldInteraction_.alpha()/(4*sqrt(pi)*pi*ewaldInteraction_.epsilon());
   }

   /*
   * Save internal state to an archive.
   */
   void MdEwaldPotential::save(Serializable::OArchive &ar)
   {
      ewaldInteraction_.save(ar);
      ar << kSpaceCutoff_;
   }

   /*
   * Get cutfoff wavenumber for long range interaction.
   */
   inline int MdEwaldPotential::nWave() const
   {  return waves_.size(); }

   /*
   * Generate waves using kCutoff; allocate memories for associated variables.
   * Generate vectors for image cell in real sapce using rSpaceCutoff. 
   * Comments:
   *   Only half of waves are stored.
   *   The first index is always non-negative.
   *
   */
   void MdEwaldPotential::makeWaves()
   {
      Vector    b0, b1, b2;    // Recprocal basis vectors.
      Vector    q0, q1, q2, q; // Partial and complete wavevectors.
      Vector    kv;            // Wavevector (as real vector)
      // double    prefactor(-0.25/ewaldInteraction_.alpha()/ewaldInteraction_.alpha());
      // double    kCutoffSq(ewaldInteraction_.kSpaceCutoffSq());
      double    ksq;
      IntVector maxK, k;       // Max and running wave indices.
      int       mink1, mink2;  // Minimum k-indices
      int       j;

      b0 = boundaryPtr_->reciprocalBasisVector(0);
      b1 = boundaryPtr_->reciprocalBasisVector(1);
      b2 = boundaryPtr_->reciprocalBasisVector(2);

      // Get max wave indices and reserve arrays
      double pi2 = 2.0*Constants::Pi;
      for (j=0; j < Dimension; ++j) {
         maxK[j] = 
             ceil(kSpaceCutoff_*boundaryPtr_->bravaisBasisVector(j).abs()/pi2);
         UTIL_CHECK(maxK[j] > 0);
      }

      if (waves_.capacity() == 0) {
         int capacity; 
         capacity = ((2*maxK[0] + 1)*(2*maxK[1] + 1)*(2*maxK[2] + 1) - 1)/2;
         UTIL_CHECK(capacity > 0);
         waves_.reserve(capacity);
         ksq_.reserve(capacity);
         g_.reserve(capacity);
         rho_.reserve(capacity);
         fexp0_.reserve(maxK[0] + 1);
         fexp1_.reserve(2*maxK[1] + 1);
         fexp2_.reserve(2*maxK[2] + 1);
      } else {
         waves_.clear();
         ksq_.clear();
         g_.clear();
         rho_.clear();
         fexp0_.clear();
         fexp1_.clear();
         fexp2_.clear();
      }

      double kCutoffSq = ewaldInteraction_.kSpaceCutoffSq();
      double pi = Constants::Pi;
      double alpha = ewaldInteraction_.alpha();
      double epsilon = ewaldInteraction_.epsilon();
      double prefactor = -0.25/(alpha*alpha);
      selfPrefactor_ = alpha/(4*sqrt(pi)*pi*epsilon);

      // Accumulate waves, and wave-related properties.
      base0_ = 0;
      upper0_ = -maxK[0];

      base1_ = maxK[1];
      upper1_ = -base1_;

      base2_ = maxK[2];
      upper2_ = -base2_;

      q0.multiply(b0, -1);
      for (k[0] = 0; k[0] <= maxK[0]; ++k[0]) { 

         // Note: First index always non-negative.
         q0 += b0;

         mink1 = (k[0] == 0 ? 0 : -maxK[1]);
         q1.multiply(b1, mink1 - 1);
         q1 += q0;
         for (k[1] = mink1; k[1] <= maxK[1]; ++k[1]) {
            q1 += b1;

            mink2 = (k[0] == 0 && k[1] == 0 ? 1 : -maxK[2]);
            q.multiply(b2, mink2 - 1);
            q += q1;

            for (k[2] = mink2; k[2] <= maxK[2]; ++k[2]) {
               q += b2;

               ksq = double(q.square());
               if (ksq <= kSpaceCutoffSq_) {

                  if (k[0] > upper0_) upper0_ = k[0];

                  if (k[1] < base1_ ) base1_  = k[1];
                  if (k[1] > upper1_) upper1_ = k[1];

                  if (k[2] < base2_ ) base2_  = k[2];
                  if (k[2] > upper2_) upper2_ = k[2];

                  waves_.append(k);
                  ksq_.append(ksq);
                  g_.append(exp(prefactor*ksq)/ksq);
               }

            } // for k[2]
         } // for k[1]
      } // for k[0]

      // Resize work arrays
      UTIL_CHECK(waves_.size() > 0);
      UTIL_CHECK(upper0_ - base0_ + 1 > 0);
      UTIL_CHECK(upper1_ - base1_ + 1 > 0);
      UTIL_CHECK(upper2_ - base2_ + 1 > 0);
      rho_.resize(waves_.size());
      fexp0_.resize(upper0_ - base0_ + 1);
      fexp1_.resize(upper1_ - base1_ + 1);
      fexp2_.resize(upper2_ - base2_ + 1);

      // Mark waves as updated
      hasWaves_ = true;
   }

   /*
   * Calculate fourier modes of charge density.
   */
   void MdEwaldPotential::computeKSpaceCharge()
   {
      // Compute waves if necessary
      if (!hasWaves()) {
         makeWaves();
      }

      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      Vector  rg;     // scaled atom position vector
      IntVector q;    // wave indices
      Vector  qv;     // wave vector
      double  dotqr;  // dot product between q and r
      double  x, y;   // real and imaginary parts of phasor
      double  charge; // atom charge
      double  TwoPi;  // 2.0*pi
      int  i, j;      // array index

      TwoPi = 2.0*Constants::Pi;

      // Clear rho for all waves
      for (i = 0; i < rho_.size(); ++i) {
         rho_[i] = DCMPLX(0.0, 0.0);
      }

      // Loop over species, molecules atoms
      int  nSpecies = simulationPtr_->nSpecies();
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         systemPtr_->begin(iSpecies, molIter); 
         for ( ; molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               boundaryPtr_->transformCartToGen(atomIter->position(), rg);
               charge = (*atomTypesPtr_)[atomIter->typeId()].charge();

               // Loop over waves
               for (i = 0; i < waves_.size(); ++i) {
                  q = waves_[i];
                  for (j = 0; j < Dimension; ++j) {
                     qv[j] = double(q[j]);
                  }
                  dotqr = rg.dot(qv)*TwoPi;
                  x = charge*cos(dotqr);
                  y = charge*sin(dotqr);
                  rho_[i] += DCMPLX(x, y);
               }

            } // For atoms.
         } // For molecules.
      } // For species.

   }

   /*
   * Add k-space Coulomb forces for all atoms.
   */
   void MdEwaldPotential::addForces()
   {
      System::MoleculeIterator molIter, imolIter, jmolIter;
      Molecule::AtomIterator atomIter, iatomIter, jatomIter;
      Vector  b[Dimension];   // array of reciprocal basis vectors
      Vector  rg;             // scaled atom position vector
      Vector  fg;             // generalized force on atom
      Vector  df;             // force contribution
      IntVector q;            // wave index
      double  charge;         // atom charge
      double  x, y;           // Real and imaginary parts of phasor
      double  TwoPi;          // 2.0*pi
      double  forcePrefactor; // 
 
      DCMPLX  TwoPiIm;       // 2.0*pi*I
      DCMPLX  de, expfactor;
      double  EPS(1.0E-10);  // Tiny number to check if is charge
      double  epsilon = ewaldInteraction_.epsilon();
      double  volume = boundaryPtr_->volume();
      int  nSpecies(simulationPtr_->nSpecies());
      int  type;
      int  i, j;

      // Constants
      TwoPi   = 2.0*Constants::Pi;
      TwoPiIm = TwoPi * Constants::Im;
      // forcePrefactor = -2.0 /(ewaldInteraction_.epsilon()*boundaryPtr_->volume());
      forcePrefactor = -2.0 /(epsilon*volume);

      // Compute Fourier components of charge density.
      computeKSpaceCharge();
    
      // Store reciprocal lattice vectors 
      for (j = 0; j < Dimension; ++j) {
         b[j] = boundaryPtr_->reciprocalBasisVector(j);
      }

      // Loop over species, molecules atoms
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         systemPtr_->begin(iSpecies, molIter); 
         for ( ; molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               type = atomIter->typeId();
               charge = (*atomTypesPtr_)[type].charge();

               if (fabs(charge) > EPS) {
                  boundaryPtr_->transformCartToGen(atomIter->position(), rg);

                  // Tabulate the exponential factors.
                  fexp0_[0] = exp(TwoPiIm * rg[0] * double(base0_));
                  de = exp(TwoPiIm * rg[0]);
                  for (i = 1; i < fexp0_.size(); ++i) {
                     fexp0_[i] = fexp0_[i-1] * de;
                  }

                  fexp1_[0] = exp(TwoPiIm * rg[1] * double(base1_));
                  de = exp(TwoPiIm * rg[1]);
                  for (i = 1; i < fexp1_.size(); ++i) {
                     fexp1_[i] = fexp1_[i-1] * de;
                  }

                  fexp2_[0] = exp(TwoPiIm * rg[2] * double(base2_));
                  de = exp(TwoPiIm * rg[2]);
                  for (i = 1; i < fexp2_.size(); ++i) {
                     fexp2_[i] = fexp2_[i-1] * de;
                  }

                  // Accumulating forces.
                  fg.zero();
                  for (i = 0; i < waves_.size(); ++i) {
                     q = waves_[i];
                     for (j = 0; j < Dimension; ++j) {
                        df[j] = double(q[j]);
                     }
                     expfactor = fexp0_[q[0]-base0_] * fexp1_[q[1]-base1_] * fexp2_[q[2]-base2_];
                     x = expfactor.real();
                     y = expfactor.imag();
                     df *= g_[i]*( x*rho_[i].imag() - y*rho_[i].real() );
                     fg += df;
                  }
                  fg *= charge*forcePrefactor;
 
                  // Transform to Cartesian coordinates
                  for (j = 0; j < Dimension; ++j) {
                     df.multiply(b[j], fg[j]);
                     atomIter->force() += df;
                  }
               }

            } // for atoms
         } // for molecules
      } // for species
   }

   /*
   * Calculate the k-space contribution to the Coulomb energy.
   */
   void MdEwaldPotential::computeEnergy()
   {
      System::MoleculeIterator imolIter, jmolIter;
      Molecule::AtomIterator iatomIter, jatomIter;
      double kPart = 0.0;
      double x, y,rhoSq;
      int nSpecies(simulationPtr_->nSpecies());


      // Compute Fourier components of charge density.
      computeKSpaceCharge();

      for (int i = 0; i < waves_.size(); ++i) {
         x = rho_[i].real();
         y = rho_[i].imag();
         rhoSq = x*x + y*y;
         kPart += ewaldInteraction_.kSpacePotential(rhoSq, g_[i]);
      }
      kPart *= 0.5 / (ewaldInteraction_.epsilon()*boundaryPtr_->volume());

      // calculate selfnergy part in ewald summation.
      double selfenergy(0.0); //store the self part energy
      double icharge;
      int iAtomType;
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         systemPtr_->begin(iSpecies, imolIter);
         for ( ; imolIter.notEnd(); ++imolIter) {
            for (imolIter->begin(iatomIter); iatomIter.notEnd(); ++iatomIter) {
               iAtomType = iatomIter->typeId();
               icharge = (*atomTypesPtr_)[iAtomType].charge();
               selfenergy += icharge * icharge;
            }
         }
      } 
      selfenergy *= selfPrefactor_;

      // Correct for conjugate wave contribution in k-part.
      kSpaceEnergy_.set(2.0 * kPart - selfenergy);
   }

   /*
   * Compute the k contribution to stress.
   * computeStress() for coulomb part in MdSystem.cpp is off.
   */
   void MdEwaldPotential::computeStress()
   {
      // Compute Fourier components of charge density.
      computeKSpaceCharge();
    
      int i;
      double x, y; //real and image part of rho[i].
      Tensor K,stressTensor;// temp stress tensor.
      IntVector q;//vector indices.
      Vector qv,q0,q1,q2; //reciprocalVector.
      Vector b0, b1, b2; // reciprocalBasisVector.
      b0 = boundaryPtr_->reciprocalBasisVector(0);
      b1 = boundaryPtr_->reciprocalBasisVector(1);
      b2 = boundaryPtr_->reciprocalBasisVector(2);

      //add a new method in Tensor.h for acquiring unit matrix.
      qv.zero();
      stressTensor.zero();

      double alpha = ewaldInteraction_.alpha();
      double ca = 0.25/(alpha*alpha);
      for (i = 0; i < waves_.size(); ++i) {
         q = waves_[i];
         q0.multiply(b0, q[0]);
         q1.multiply(b1, q[1]);
         q2.multiply(b2, q[2]); 
         qv += q0;
         qv += q1;
         qv += q2;

         K.dyad(qv,qv);
         K *=  -2.0 * (ca + (1.0 / ksq_[i]));
         K.add(unitTensor_, K);
         x = rho_[i].real();
         y = rho_[i].imag();
         K *= g_[i]*( x*x + y*y);
         stressTensor += K;  
      }

      double epsilon = ewaldInteraction_.epsilon();
      double volume = boundaryPtr_->volume(); 
      stressTensor /= (2.0 * epsilon * volume);   

      // Correct for conjugate wave contribution in k-part then set.
      stressTensor *= 2.0; 

      kSpaceStress_.set(stressTensor);
   }

} 
#endif
