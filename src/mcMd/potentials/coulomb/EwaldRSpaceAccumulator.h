#ifndef MD_EWALD_RSPACE_ACCUMULATOR_h
#define MD_EWALD_RSPACE_ACCUMULATOR_h 


#include <util/space/Tensor.h>
#include <util/misc/Setable.h>

namespace McMd
{

  using namespace Util;

  class EwaldRSpaceAccumulator 
  {

     template <typename T>
     friend class MdEwaldPairPotentialImpl;

     public:

     /*
     * Constructor.
     */
     EwaldRSpaceAccumulator()
        :rSpaceEnergy_(),
         rSpaceStress_()
     {
        rSpaceEnergy_.unset();
        rSpaceStress_.unset();
     }

     bool isSetEnergy() const
     {  return rSpaceEnergy_.isSet(); }

     bool isSetStress() const
     { return rSpaceStress_.isSet(); }

     double rSpaceEnergy() const
     { return rSpaceEnergy_.value(); }

     Tensor rSpaceStress() const
     { return rSpaceStress_.value(); }
     
     //Tensor rSpacePressure() const;

     protected:
     Setable<double> rSpaceEnergy_;
     Setable<Tensor> rSpaceStress_;
                  
  };
}
#endif
