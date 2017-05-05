#ifndef MD_EWALD_POTENTIAL_CPP
#define MD_EWALD_POTENTIAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "MdPMEPotential.h" 
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <util/space/Vector.h>
#include <util/space/Tensor.h>
#include <util/math/Constants.h>
#include <stdlib.h>

#include <util/boundary/Boundary.h>
#include <cmath>
#include <util/containers/Array.h>
#include <util/containers/DArray.h>
#include <mcMd/chemistry/AtomType.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   MdPMEPotential::MdPMEPotential(System& system)
    : MdCoulombPotential(),
      rSpaceAccumulator_(),
      ewaldInteraction_(),
      simulationPtr_(&system.simulation()),
      systemPtr_(&system),
      boundaryPtr_(&system.boundary()),
      atomTypesPtr_(&system.simulation().atomTypes()),
      gridSize_(),
      Qgrid_(),
      Qhatgrid_(),
      BCgrid_(),
      ikop_(),
      xfield_(),
      yfield_(),
      zfield_(),
      order_(5),
      BCikinitialized_(false)
   {
      // initialize unit tensor.
      const double unitMatrix_[3][3] = { {1,0,0}, {0,1,0}, {0,0,1}};
      Tensor unitTensor_(unitMatrix_);
      // Note: Don't setClassName - using "CoulombPotential" base class name
   }

   /*
   * Destructor (does nothing)
   */
   MdPMEPotential::~MdPMEPotential()
   {
      fftw_destroy_plan(forward_plan);
      fftw_destroy_plan(xfield_backward_plan);
      fftw_destroy_plan(yfield_backward_plan);
      fftw_destroy_plan(zfield_backward_plan);
   }

   /*
   * Read parameters and initialize.
   */
   void MdPMEPotential::readParameters(std::istream& in)
   {
      ewaldInteraction_.readParameters(in);
      read<int>(in, "xgridSize", gridSize_[0]);
      read<int>(in, "ygridSize", gridSize_[1]);
      read<int>(in, "zgridSize", gridSize_[2]);

      // allocate memory for grid
      Qgrid_.allocate(gridSize_);
      Qhatgrid_.allocate(gridSize_);
      BCgrid_.allocate(gridSize_);
      ikop_.allocate(std::max(std::max(gridSize_[1],gridSize_[2]),gridSize_[0]));

      xfield_.allocate(gridSize_);
      yfield_.allocate(gridSize_);
      zfield_.allocate(gridSize_);

      // initialize fft plan for Q grid.
      fftw_complex* inf  = reinterpret_cast<fftw_complex*>(Qgrid_.data());
      fftw_complex* outf = reinterpret_cast<fftw_complex*>(Qhatgrid_.data());
      forward_plan = fftw_plan_dft_3d(gridSize_[0],gridSize_[1],gridSize_[2],
                                      inf, outf, FFTW_FORWARD,FFTW_MEASURE);

      // initialize fft plan for electric field grid.
      fftw_complex* inxf  = reinterpret_cast<fftw_complex*>(xfield_.data());
      fftw_complex* outxf = reinterpret_cast<fftw_complex*>(xfield_.data());
      xfield_backward_plan = fftw_plan_dft_3d(gridSize_[0],gridSize_[1],gridSize_[2],
                                      inxf, outxf,FFTW_BACKWARD, FFTW_MEASURE);
 
      fftw_complex* inyf  = reinterpret_cast<fftw_complex*>(yfield_.data());
      fftw_complex* outyf = reinterpret_cast<fftw_complex*>(yfield_.data());
      yfield_backward_plan = fftw_plan_dft_3d(gridSize_[0],gridSize_[1],gridSize_[2],
                                      inyf, outyf,FFTW_BACKWARD, FFTW_MEASURE);
    
      fftw_complex* inzf  = reinterpret_cast<fftw_complex*>(zfield_.data());
      fftw_complex* outzf = reinterpret_cast<fftw_complex*>(zfield_.data());
      zfield_backward_plan = fftw_plan_dft_3d(gridSize_[0],gridSize_[1],gridSize_[2],
                                      inzf, outzf,FFTW_BACKWARD, FFTW_MEASURE);
 
      //Calculate prefactors frequently used in this class.
      double pi = Constants::Pi;
      selfPrefactor_ = ewaldInteraction_.alpha()/(4*sqrt(pi)*pi*ewaldInteraction_.epsilon());
   }

   /*
   * virtual function. place holder.
   */
   inline int MdPMEPotential::nWave() const
   {  return 0; }

   /*
   * virtual function. place holder. 
   */
   void MdPMEPotential::makeWaves()
   {; }

   /*
   * set elements of grid to all zero.
   */
   template<class T>
   void MdPMEPotential::initializeGrid(GridArray<T>& grid) 
   {
      IntVector temp;

      for (int i = 0; i < grid.dimension(0); ++i) {
         temp[0] = i;
         for (int j = 0; j < grid.dimension(1); ++j) {
            temp[1] = j;
            for (int k = 0; k < grid.dimension(2); ++k) {
               temp[2] = k;
               grid(temp) = 0 ;
            }
         }
      }
   }

   /*
   * BCgrid_ --- BC(g1,g2,g3).
   */
   void MdPMEPotential::influence_function()
   {
      double m0, m1, m2; //temp parameter for C grid.
      double msqr;
      Vector b0, b1, b2;
      double b,c;
      IntVector gridPoint; //grid point coordination.
      double pi(Constants::Pi);
      double pi2(pi*pi);
      double alpha(ewaldInteraction_.alpha());
      double alpha2(alpha*alpha);

      initializeGrid(BCgrid_);

      b0 = boundaryPtr_->reciprocalBasisVector(0) ;
      b1 = boundaryPtr_->reciprocalBasisVector(1) ;
      b2 = boundaryPtr_->reciprocalBasisVector(2) ;


      // Loop over i, j, k grid point.
      for ( double i = 0.0; i < gridSize_[0]; ++i) {
         gridPoint[0] = i;
         m0 = (i <= gridSize_[0] / 2.0) ? i : i - gridSize_[0];
         for ( double j = 0.0; j < gridSize_[1]; ++j) {
            gridPoint[1] = j;
            m1 = (j <= gridSize_[1] / 2.0) ? j : j - gridSize_[1];
            for ( double k = 0.0; k < gridSize_[2]; ++k) {
               gridPoint[2] = k;
               m2 = (k <= gridSize_[2] / 2.0) ? k : k - gridSize_[2];

               b0 = boundaryPtr_->reciprocalBasisVector(0) ;
               b1 = boundaryPtr_->reciprocalBasisVector(1) ;
               b2 = boundaryPtr_->reciprocalBasisVector(2) ;

               b0 *= m0/(2.0*pi);
               b1 *= m1/(2.0*pi);
               b2 *= m2/(2.0*pi);

               b0 += b1;
               b0 += b2;

               msqr   = b0.dot(b0);

               b = bfactor(i,0) * bfactor(j,1) *bfactor(k,2);
               c = 1.0 / (4*pi2* ewaldInteraction_.epsilon() * boundaryPtr_->volume() * msqr)
                 * exp( -pi2*msqr/alpha2);

               BCgrid_(gridPoint) = b * c;
            }
         }
      }
      BCgrid_[0] = 0;
   }
 
   /*
   * component of B grid
   */
   double MdPMEPotential::bfactor(double m, int dim)
   {
      double pi(Constants::Pi);
      double gridSize(gridSize_[dim]);
      DCMPLX I(0.0,1.0);

      DCMPLX denom(0.0, 0.0);
      for (double k = 0.0; k <= order_ - 2.0; k++) {
         denom += basisSpline(k + 1.0)*exp(2.0 * pi * I * m * k / gridSize) ;
      }
      return std::norm(exp(2.0 * pi * I * (order_ -1.0) * m/gridSize) / denom);
   }
 
   /*
   * Qgrid_ --- Q(g1,g2,g3).
   */
   void MdPMEPotential::spreadCharge()
   {
      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      double  charge;
      Vector  gpos; //general coordination of atom
      int ximg, yimg, zimg; // grid point coordination
      double xdistance, ydistance, zdistance;//distance between atom and grid point.
      int xknot, yknot, zknot;
      IntVector knot;
      IntVector floorGridIdx;

      initializeGrid(Qgrid_);

      // Loop over atoms.
      int  nSpecies = simulationPtr_->nSpecies();
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         systemPtr_->begin(iSpecies, molIter);
         for ( ; molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               charge = (*atomTypesPtr_)[atomIter->typeId()].charge();
               boundaryPtr_->transformCartToGen(atomIter->position(), gpos);

               // Find the floor grid point.
               floorGridIdx[0]=floor(gpos[0]*gridSize_[0]);
               floorGridIdx[1]=floor(gpos[1]*gridSize_[1]);
               floorGridIdx[2]=floor(gpos[2]*gridSize_[2]);

               // haven't incoorporate reciprocal vector * lattice vector in other kind lattice.
               // need to modify the expression of  distance.
               for (int x = 0 ; x < order_ ; ++x) {
                  ximg = floorGridIdx[0] + x - (order_ - 1);
                  xdistance = (gpos[0]*gridSize_[0]-ximg) ;
                  xknot = ximg < 0 ? ximg + gridSize_[0] : ximg;
                  knot[0] = xknot;

                  for (int y = 0 ; y < order_ ; ++y) {
                     yimg = floorGridIdx[1] + y - (order_ - 1);
                     ydistance = (gpos[1]*gridSize_[1]-yimg) ;
                     yknot =  yimg < 0 ? yimg + gridSize_[1] : yimg;
                     knot[1] = yknot;

                     for (int z = 0; z < order_; ++z) {
                        zimg = floorGridIdx[2] + z - (order_ - 1);
                        zdistance = (gpos[2]*gridSize_[2]-zimg) ;
                        zknot =  zimg < 0 ? zimg + gridSize_[2] : zimg;
                        knot[2] = zknot;

                        Qgrid_(knot) += charge 
                                      * basisSpline(xdistance)
                                      * basisSpline(ydistance)
                                      * basisSpline(zdistance);
                     }//loop z
                  }//loop y
               }//loop x
            }//loop atom
         }//loop molecule
      }//loop over species
   } 
        
   /*
   * basisSpline function.
   */
   double MdPMEPotential::basisSpline(double x)
   {
      if (0.0 < x && x < 1.0)
         return 1.0 / 24.0 * x * x * x *x;
      else if (1.0 <= x && x < 2.0)
         return -5.0/24.0 + (5.0/6.0 + (-5.0/4.0 + (5.0 / 6.0 -1.0/6.0*x) *x) *x) *x;
      else if (2.0 <= x && x < 3.0)
         return 155.0/24.0 + (-25.0/2.0 + (35.0/4.0 + (-5.0/2.0 +1.0/4.0*x)*x)*x)*x;
      else if (3.0 <= x && x < 4.0)
         return -655.0/24.0 + (65.0/2.0 + ( -55.0/4.0 + (5.0/2.0 - 1.0/6.0*x)*x)*x)*x;
      else if (4.0 <= x && x < 5.0)
         return 625.0/24.0 + (-125.0/6.0 + (25.0/4.0 + (-5.0/6.0 + 1/24.0 *x)*x)*x)*x;
      else
         return 0.0;
   }
   
   /*
   * Calculate the k-space contribution to the Coulomb energy.
   */
   void MdPMEPotential::computeEnergy()
   {
      // unset energy accumulator.
      //unsetEnergy();

      System::MoleculeIterator imolIter, jmolIter;
      Molecule::AtomIterator iatomIter, jatomIter;
      double fft_k_energy = 0.0;
      int nSpecies(simulationPtr_->nSpecies());
      IntVector pos;

      spreadCharge();
 
      if(!BCikinitialized_){
         influence_function();
         ik_differential_operator();
         BCikinitialized_=true;
      }

      initializeGrid(Qhatgrid_);
      fftw_execute(forward_plan);

      for (int i = 0; i < gridSize_[0]; ++i) {
         for (int j = 0; j < gridSize_[1]; ++j) {
            for (int k = 0; k < gridSize_[2]; ++k) {
               pos[0] = i;
               pos[1] = j;
               pos[2] = k;

               fft_k_energy += BCgrid_(pos) * std::norm(Qhatgrid_(pos));
            }
         }
      }

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

      kSpaceEnergy_.set(0.5*fft_k_energy - selfenergy);
  }

   // n-level rather than k level ie. missing the prefactor 2.0 * Pi * I / L. 
   void MdPMEPotential::ik_differential_operator()
   {
      for (int i = 0 ; i < 3 ; ++i){
         ikop_[0][i] = 0.0;
         ikop_[gridSize_[i]/2][i] = 0.0;
         for (int j = 1; j < ikop_.capacity()/2; ++j) {
            ikop_[j][i] = j;
            ikop_[gridSize_[i] - j][i] = -j;
         }
      }
   }

   /*
   * Add k-space Coulomb forces for all atoms. standard Ewald Summation.
   */
   void MdPMEPotential::addForces()
   {
      System::MoleculeIterator molIter, imolIter, jmolIter;
      Molecule::AtomIterator atomIter, iatomIter, jatomIter;
      Vector fatom(0.0);
      Vector floorGridIdx, gpos;
      IntVector m;
      IntVector knot;
      double xdistance, ydistance, zdistance; //b-spline 
      int ximg, yimg, zimg;
      int xknot, yknot, zknot;
      double  charge;         // atom charge
      double  TwoPi;          // 2.0*pi
      double  forcePrefactor; // 
      DCMPLX  TwoPiIm;       // 2.0*pi*I
      double  EPS(1.0E-10);  // Tiny number to check if is charge
      int  nSpecies(simulationPtr_->nSpecies());
      int  type;
      int  pos; //rank in GridArray


      // Constants
      TwoPi   = 2.0*Constants::Pi;
      TwoPiIm = TwoPi * Constants::Im;
      forcePrefactor = -2.0 /(ewaldInteraction_.epsilon()*boundaryPtr_->volume());

      spreadCharge();
 
      if(!BCikinitialized_){
         influence_function();
         ik_differential_operator();
         BCikinitialized_=true;
      }
 
      initializeGrid(Qhatgrid_);
 
      fftw_execute(forward_plan);
      for (int i = 0 ; i < gridSize_[0] ; ++i) {
         for (int j = 0 ; j < gridSize_[1] ; ++j) {
            for (int k = 0 ; k < gridSize_[2] ; ++k) {

               pos = i * gridSize_[1]*gridSize_[2] + j * gridSize_[2] + k;

               xfield_[pos] = TwoPiIm / boundaryPtr_->length(0) 
                            * double(ikop_[i][0]) * Qhatgrid_[pos] * BCgrid_[pos];
               yfield_[pos] = TwoPiIm / boundaryPtr_->length(1) 
                            * double(ikop_[j][1]) * Qhatgrid_[pos] * BCgrid_[pos];
               zfield_[pos] = TwoPiIm / boundaryPtr_->length(2) 
                            * double(ikop_[k][2]) * Qhatgrid_[pos] * BCgrid_[pos];
            }//loop over z
         }//loop over y
      }//loop over x

      fftw_execute(xfield_backward_plan);
      fftw_execute(yfield_backward_plan);
      fftw_execute(zfield_backward_plan);
      // Loop over species, molecules atoms
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         systemPtr_->begin(iSpecies, molIter); 
         for ( ; molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               type = atomIter->typeId();
               charge = (*atomTypesPtr_)[type].charge();

               if( fabs(charge) > EPS) {
                  boundaryPtr_->transformCartToGen(atomIter->position(), gpos);

                  // Find the floor grid point.
                  floorGridIdx[0]=floor(gpos[0]*gridSize_[0]);
                  floorGridIdx[1]=floor(gpos[1]*gridSize_[1]);
                  floorGridIdx[2]=floor(gpos[2]*gridSize_[2]);
    
                  fatom.zero();
                  for (int x = 0 ; x < order_ ; ++x) {
                     ximg = floorGridIdx[0] + x - (order_ - 1);
                     xdistance = (gpos[0]*gridSize_[0]-ximg) ;
                     xknot = ximg < 0 ? ximg + gridSize_[0] : ximg;
                     knot[0] = xknot;
   
                     for (int y = 0 ; y < order_ ; ++y) {
                        yimg = floorGridIdx[1] + y - (order_ - 1);
                        ydistance = (gpos[1]*gridSize_[1]-yimg) ;
                        yknot =  yimg < 0 ? yimg + gridSize_[1] : yimg;
                        knot[1] = yknot;
   
                        for (int z = 0; z < order_; ++z) {
                           zimg = floorGridIdx[2] + z - (order_ - 1);
                           zdistance = (gpos[2]*gridSize_[2]-zimg) ;
                           zknot =  zimg < 0 ? zimg + gridSize_[2] : zimg;
                           knot[2] = zknot;
   
                           fatom[0] +=basisSpline(xdistance)
                                     *basisSpline(ydistance)
                                     *basisSpline(zdistance)
                                     *std::real(xfield_(knot));
                           fatom[1] +=basisSpline(xdistance)
                                     *basisSpline(ydistance)
                                     *basisSpline(zdistance)
                                     *std::real(yfield_(knot));
                           fatom[2] +=basisSpline(xdistance)
                                     *basisSpline(ydistance)
                                     *basisSpline(zdistance)
                                     *std::real(zfield_(knot));
                        }
                     }
                  }
                  fatom *= -1.0*charge;
                  atomIter->force() += fatom;
               }// if charged
            }// loop over atom
         }// loop over mol 
      }// loop over species
   }

   /*
   * place holder 
   */
   void MdPMEPotential::computeStress()
   { ;}
} 
#endif
