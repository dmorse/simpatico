#ifndef HARMONIC_ANGLE_TEST_H
#define HARMONIC_ANGLE_TEST_H

#include <simp/interaction/angle/HarmonicAngle.h>
#include "AngleTestTemplate.h"

class HarmonicAngleTest : public AngleTestTemplate<HarmonicAngle>
{
public:

   using AngleTestTemplate<HarmonicAngle>::readParamFile;

   void setUp()
   {
      eps_ = 1.0E-5;
      setNAngleType(1);  
      readParamFile("in/HarmonicAngle");  
   } 

   void testSetUp()
   {  printMethod(TEST_FUNC); } 

   void testForce()
   {
      printMethod(TEST_FUNC);
      b1_ = Vector( 0.2, 0.90,  0.3);
      b2_ = Vector(-0.1, 0.85, -0.4);
      type_ = 0;
      forceTest();
   }

};

TEST_BEGIN(HarmonicAngleTest)
TEST_ADD(HarmonicAngleTest, testSetUp)
TEST_ADD(HarmonicAngleTest, testForce)
TEST_END(HarmonicAngleTest)

#endif
