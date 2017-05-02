#ifndef COSINE_ANGLE_TEST_H
#define COSINE_ANGLE_TEST_H

#include <simp/interaction/angle/CosineAngle.h>
#include "AngleTestTemplate.h"

class CosineAngleTest : public AngleTestTemplate<CosineAngle>
{
public:

   #if 0
   CosineAngleTest() :
      AngleTestTemplate<CosineAngle>("in/CosineAngle")
   {}
   #endif

   using AngleTestTemplate<CosineAngle>::readParamFile;

   void setUp()
   {
      eps_ = 1.0E-5;
      setNAngleType(1);  
      readParamFile("in/CosineAngle");  
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

TEST_BEGIN(CosineAngleTest)
TEST_ADD(CosineAngleTest, testSetUp)
TEST_ADD(CosineAngleTest, testForce)
TEST_END(CosineAngleTest)

#endif
