/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "MdSpmePotential.h" 
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
   MdSpmePotential::MdSpmePotential(System& system)
    : MdCoulombPotential(),
      ewaldInteraction_(),
      simulationPtr_(&system.simulation()),
      systemPtr_(&system),
      boundaryPtr_(&system.boundary()),
      atomTypesPtr_(&system.simulation().atomTypes()),
      gridDimensions_(),
      Qgrid_(),
      Qhatgrid_(),
      BCgrid_(),
      ikop_(),
      xfield_(),
      yfield_(),
      zfield_(),
      order_(5)
   {
      //initialize unit tensor.
      const double unitMatrix_[3][3] = { {1,0,0}, {0,1,0}, {0,0,1}};
      Tensor unitTensor_(unitMatrix_);
      // Note: Don't setClassName - using "CoulombPotential" base class name
   }

   /*
   * Destructor 
   */
   MdSpmePotential::~MdSpmePotential()
   {
      fftw_destroy_plan(forward_plan);
      fftw_destroy_plan(xfield_backward_plan);
      fftw_destroy_plan(yfield_backward_plan);
      fftw_destroy_plan(zfield_backward_plan);
   }

   /*
   * Read parameters and initialize.
   */
   void MdSpmePotential::readParameters(std::istream& in)
   {
      // Read EwaldInteraction block containing parameters
      bool nextIndent = false;
      addParamComposite(ewaldInteraction_, nextIndent);
      ewaldInteraction_.readParameters(in);
      read<IntVector>(in, "gridDimensions",gridDimensions_);
   }

   /*
   * Load internal state from an archive.
   */
   void MdSpmePotential::loadParameters(Serializable::IArchive &ar)
   {
      bool nextIndent = false;
      addParamComposite(ewaldInteraction_, nextIndent);
      ewaldInteraction_.loadParameters(ar);
      loadParameter<IntVector>(ar, "gridDimensions", gridDimensions_);
   }

   /*
   * Save internal state to an archive.
   */
   void MdSpmePotential::save(Serializable::OArchive &ar)
   {
      ewaldInteraction_.save(ar);
      ar << gridDimensions_;
   }

   /*
   * Set a parameter value, identified by a string.
   */
   void MdSpmePotential::set(std::string name, double value)
   {
      ewaldInteraction_.set(name, value); 
      unsetWaves();
   }

   /*
   * Get a parameter value, identified by a string.
   */
   double MdSpmePotential::get(std::string name) const
   {
      double value;
      value = ewaldInteraction_.get(name); 
      return value;
   }

   /*
   * place holder 
   */
   inline int MdSpmePotential::nWave() const
   {  return gridDimensions_[0]*gridDimensions_[1]*gridDimensions_[2]; }

   /*
   * Precompute waves and influence function.
   */
   void MdSpmePotential::makeWaves()
   { 
      if (BCgrid_.size() == 0) {

         // Allocate memory for grid
         Qgrid_.allocate(gridDimensions_);
         Qhatgrid_.allocate(gridDimensions_);
         BCgrid_.allocate(gridDimensions_);
         int maxDim = std::max(gridDimensions_[1],gridDimensions_[2]);
         maxDim = std::max(maxDim, gridDimensions_[0]);
         ikop_.allocate(maxDim);
         xfield_.allocate(gridDimensions_);
         yfield_.allocate(gridDimensions_);
         zfield_.allocate(gridDimensions_);

         // Initialize fft plan for Q grid
         fftw_complex* inf;
         fftw_complex* outf;
         inf  = reinterpret_cast<fftw_complex*>(Qgrid_.data());
         outf = reinterpret_cast<fftw_complex*>(Qhatgrid_.data());
         forward_plan = fftw_plan_dft_3d(gridDimensions_[0],
                                         gridDimensions_[1],
                                         gridDimensions_[2],
                                         inf, outf, 
                                         FFTW_FORWARD,FFTW_MEASURE);

         // Initialize fft plans for electric field component grids
         fftw_complex* inxf; 
         fftw_complex* outxf;
         inxf = reinterpret_cast<fftw_complex*>(xfield_.data());
         outxf = reinterpret_cast<fftw_complex*>(xfield_.data());
         xfield_backward_plan = 
                  fftw_plan_dft_3d(gridDimensions_[0],
                                   gridDimensions_[1],
                                   gridDimensions_[2],
                                   inxf, outxf,FFTW_BACKWARD, FFTW_MEASURE);
         fftw_complex* inyf; 
         fftw_complex* outyf;
         inyf  = reinterpret_cast<fftw_complex*>(yfield_.data());
         outyf = reinterpret_cast<fftw_complex*>(yfield_.data());
         yfield_backward_plan = 
                  fftw_plan_dft_3d(gridDimensions_[0],
                                   gridDimensions_[1],
                                   gridDimensions_[2],
                                   inyf, outyf,FFTW_BACKWARD, FFTW_MEASURE);
         fftw_complex* inzf; 
         fftw_complex* outzf;
         inzf  = reinterpret_cast<fftw_complex*>(zfield_.data());
         outzf = reinterpret_cast<fftw_complex*>(zfield_.data());
         zfield_backward_plan = 
                  fftw_plan_dft_3d(gridDimensions_[0],
                                   gridDimensions_[1],
                                   gridDimensions_[2],
                                   inzf, outzf,FFTW_BACKWARD, FFTW_MEASURE);
      }

      influence_function();
      ik_differential_operator();
  
      // Mark wave data as up to date.
      hasWaves_ = true;
   }

   /*
   * Set elements of grid to all zero.
   */
   template<class T>
   void MdSpmePotential::setGridToZero(GridArray<T>& grid) 
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
   void MdSpmePotential::influence_function()
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

      setGridToZero(BCgrid_);

      b0 = boundaryPtr_->reciprocalBasisVector(0) ;
      b1 = boundaryPtr_->reciprocalBasisVector(1) ;
      b2 = boundaryPtr_->reciprocalBasisVector(2) ;


      // Loop over i, j, k grid point.
      for ( double i = 0.0; i < gridDimensions_[0]; ++i) {
         gridPoint[0] = i;
         m0 = (i <= gridDimensions_[0] / 2.0) ? i : i - gridDimensions_[0];
         for ( double j = 0.0; j < gridDimensions_[1]; ++j) {
            gridPoint[1] = j;
            m1 = (j <= gridDimensions_[1] / 2.0) ? j : j - gridDimensions_[1];
            for ( double k = 0.0; k < gridDimensions_[2]; ++k) {
               gridPoint[2] = k;
               m2 = (k <= gridDimensions_[2] / 2.0) ? k : k - gridDimensions_[2];

               b0 = boundaryPtr_->reciprocalBasisVector(0) ;
               b1 = boundaryPtr_->reciprocalBasisVector(1) ;
               b2 = boundaryPtr_->reciprocalBasisVector(2) ;

               b0 *= m0/(2.0*pi);
               b1 *= m1/(2.0*pi);
               b2 *= m2/(2.0*pi);

               b0 += b1;
               b0 += b2;

               msqr   = b0.dot(b0);

               b = bfactor(i, 0) * bfactor(j, 1) *bfactor(k, 2);
               c = 1.0 / (4.0*pi2* ewaldInteraction_.epsilon()*msqr)
                   * exp( -pi2*msqr/alpha2);
               c /= boundaryPtr_->volume();

               BCgrid_(gridPoint) = b * c;
            }
         }
      }
      BCgrid_[0] = 0;
   }

   /*
   * component of B grid
   */
   double MdSpmePotential::bfactor(double m, int dim)
   {
      double pi(Constants::Pi);
      double gridDimensions(gridDimensions_[dim]);
      DCMPLX I(0.0,1.0);

      DCMPLX denom(0.0, 0.0);
      
      // If order of spline is odd, this interpolation result fails,
      // when 2*m = gridDimensions[dim], since 1 / 0 in this function.
      if( order_ % 2 == 1 && m == gridDimensions_[dim] / 2.0) {
         m = m-1 ;
      }

      for (double k = 0.0; k <= order_ - 2.0; k++) {
         denom += basisSpline(k + 1.0)*exp(2.0 * pi * I * m * k / gridDimensions) ;
      }
      return std::norm(exp(2.0 * pi * I * (order_ -1.0) * m/gridDimensions) / denom);
   }
 
   /*
   * Qgrid_ --- Q(g1,g2,g3).
   */
   void MdSpmePotential::spreadCharge()
   {
      if (!hasWaves()) {
         makeWaves();
      }

      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      double  charge;
      Vector  gpos; //general coordination of atom
      int ximg, yimg, zimg; // grid point coordination
      double xdistance, ydistance, zdistance;//distance between atom and grid point.
      int xknot, yknot, zknot;
      IntVector knot;
      IntVector floorGridIdx;

      setGridToZero(Qgrid_);

      // Loop over atoms.
      int  nSpecies = simulationPtr_->nSpecies();
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         systemPtr_->begin(iSpecies, molIter);
         for ( ; molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               charge = (*atomTypesPtr_)[atomIter->typeId()].charge();

               // Compute generalized position with components in [0,1]
               boundaryPtr_->transformCartToGen(atomIter->position(), gpos);
               boundaryPtr_->shiftGen(gpos);

               // Find the floor grid point.
               floorGridIdx[0]=floor(gpos[0]*gridDimensions_[0]);
               floorGridIdx[1]=floor(gpos[1]*gridDimensions_[1]);
               floorGridIdx[2]=floor(gpos[2]*gridDimensions_[2]);

               // Have not incorporated reciprocal vector * lattice vector 
               // for other lattice. Need to modify the expression of  distance.

               for (int x = 0 ; x < order_ ; ++x) {
                  ximg = floorGridIdx[0] + x - (order_ - 1);
                  xdistance = (gpos[0]*gridDimensions_[0]-ximg);
                  xknot = ximg < 0 ? ximg + gridDimensions_[0] : ximg;
                  knot[0] = xknot;

                  for (int y = 0 ; y < order_ ; ++y) {
                     yimg = floorGridIdx[1] + y - (order_ - 1);
                     ydistance = (gpos[1]*gridDimensions_[1]-yimg) ;
                     yknot =  yimg < 0 ? yimg + gridDimensions_[1] : yimg;
                     knot[1] = yknot;

                     for (int z = 0; z < order_; ++z) {
                        zimg = floorGridIdx[2] + z - (order_ - 1);
                        zdistance = (gpos[2]*gridDimensions_[2]-zimg) ;
                        zknot =  zimg < 0 ? zimg + gridDimensions_[2] : zimg;
                        knot[2] = zknot;

                        Qgrid_(knot) += charge 
                                      * basisSpline(xdistance)
                                      * basisSpline(ydistance)
                                      * basisSpline(zdistance);
                     } //loop z
                  } //loop y
               } //loop x
            } //loop atom
         } //loop molecule
      } //loop over species
   } 
 
   /*
   * basisSpline function.
   */
   double MdSpmePotential::basisSpline(double x)
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
   *  n-level rather than k level ie. missing the prefactor 2.0 * Pi * I / L. 
   */
   void MdSpmePotential::ik_differential_operator()
   {
      for (int i = 0 ; i < 3 ; ++i){
         ikop_[0][i] = 0.0;
         ikop_[gridDimensions_[i]/2][i] = 0.0;
         for (int j = 1; j < ikop_.capacity()/2; ++j) {
            ikop_[j][i] = j;
            ikop_[gridDimensions_[i] - j][i] = -j;
         }
      }
   }

   /*
   * Add k-space Coulomb forces for all atoms. standard Ewald Summation.
   */
   void MdSpmePotential::addForces()
   {
      spreadCharge();
 
      setGridToZero(Qhatgrid_);
      fftw_execute(forward_plan);

      double TwoPi   = 2.0*Constants::Pi;      // 2*Pi
      DCMPLX TwoPiIm = TwoPi * Constants::Im;  // 2*Pi*I
      int  pos;                                //rank in GridArray

      for (int i = 0 ; i < gridDimensions_[0] ; ++i) {
         for (int j = 0 ; j < gridDimensions_[1] ; ++j) {
            for (int k = 0 ; k < gridDimensions_[2] ; ++k) {

               pos = i * gridDimensions_[1]*gridDimensions_[2] + j * gridDimensions_[2] + k;

               xfield_[pos] = TwoPiIm / boundaryPtr_->length(0) 
                            * double(ikop_[i][0]) * Qhatgrid_[pos] * BCgrid_[pos];
               yfield_[pos] = TwoPiIm / boundaryPtr_->length(1) 
                            * double(ikop_[j][1]) * Qhatgrid_[pos] * BCgrid_[pos];
               zfield_[pos] = TwoPiIm / boundaryPtr_->length(2) 
                            * double(ikop_[k][2]) * Qhatgrid_[pos] * BCgrid_[pos];
            } // loop over z
         } // loop over y
      } // loop over x

      fftw_execute(xfield_backward_plan);
      fftw_execute(yfield_backward_plan);
      fftw_execute(zfield_backward_plan);

      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      Vector fatom(0.0);
      Vector floorGridIdx, gpos;
      IntVector m;
      IntVector knot;
      double xdistance, ydistance, zdistance; // b-spline 
      double  EPS = 1.0E-10;  // Tiny number to check if is charged
      double charge;  
      int type;
      int ximg, yimg, zimg;
      int xknot, yknot, zknot;

      // Loop over species, molecules atoms
      int  nSpecies(simulationPtr_->nSpecies());
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         systemPtr_->begin(iSpecies, molIter); 
         for ( ; molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               type = atomIter->typeId();
               charge = (*atomTypesPtr_)[type].charge();

               // If atom has nonzero charge 
               if( fabs(charge) > EPS) {

                  // Compute generalized position with components in [0,1]
                  boundaryPtr_->transformCartToGen(atomIter->position(), gpos);
                  boundaryPtr_->shiftGen(gpos);

                  // Find the floor grid point.
                  floorGridIdx[0]=floor(gpos[0]*gridDimensions_[0]);
                  floorGridIdx[1]=floor(gpos[1]*gridDimensions_[1]);
                  floorGridIdx[2]=floor(gpos[2]*gridDimensions_[2]);
   
                  // Compute the force on the atom by interpolation 
                  fatom.zero();
                  for (int x = 0 ; x < order_ ; ++x) {
                     ximg = floorGridIdx[0] + x - (order_ - 1);
                     xdistance = (gpos[0]*gridDimensions_[0]-ximg) ;
                     xknot = ximg < 0 ? ximg + gridDimensions_[0] : ximg;
                     knot[0] = xknot;
   
                     for (int y = 0 ; y < order_ ; ++y) {
                        yimg = floorGridIdx[1] + y - (order_ - 1);
                        ydistance = (gpos[1]*gridDimensions_[1]-yimg) ;
                        yknot =  yimg < 0 ? yimg + gridDimensions_[1] : yimg;
                        knot[1] = yknot;
   
                        for (int z = 0; z < order_; ++z) {
                           zimg = floorGridIdx[2] + z - (order_ - 1);
                           zdistance = (gpos[2]*gridDimensions_[2]-zimg) ;
                           zknot =  zimg < 0 ? zimg + gridDimensions_[2] : zimg;
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
               } // if charged
            } // loop over atom
         } // loop over mol 
      } // loop over species
   }

   /*
   * Calculate the k-space contribution to the Coulomb energy.
   */
   void MdSpmePotential::computeEnergy()
   {
      spreadCharge();
 
      setGridToZero(Qhatgrid_);
      fftw_execute(forward_plan);

      IntVector pos;
      double energy = 0.0;
      int i, j, k;
      for (i = 0; i < gridDimensions_[0]; ++i) {
         pos[0] = i;
         for (j = 0; j < gridDimensions_[1]; ++j) {
            pos[1] = j;
            for (k = 0; k < gridDimensions_[2]; ++k) {
               pos[2] = k;
               energy += BCgrid_(pos) * std::norm(Qhatgrid_(pos));
            }
         }
      }
      energy /= 2.0;

      // Calculate self-energy correction to Ewald summation.
      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      double charge;
      double selfEnergy = 0.0; 
      int nSpecies = simulationPtr_->nSpecies();
      for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
         systemPtr_->begin(iSpecies, molIter);
         for ( ; molIter.notEnd(); ++molIter) {
            molIter->begin(atomIter); 
            for ( ; atomIter.notEnd(); ++atomIter) {
               charge = (*atomTypesPtr_)[atomIter->typeId()].charge();
               selfEnergy += charge*charge;
            }
         }
      }
      double pi = Constants::Pi;
      double alpha = ewaldInteraction_.alpha();
      double epsilon = ewaldInteraction_.epsilon();
      selfEnergy *= alpha/(4.0*sqrt(pi)*pi*epsilon);
      energy -= selfEnergy;

      kSpaceEnergy_.set(energy);
   }

   /*
   * place holder 
   */
   void MdSpmePotential::computeStress()
   { ;}

} 
