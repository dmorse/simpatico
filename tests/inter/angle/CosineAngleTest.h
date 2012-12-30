#ifndef COSINE_ANGLE_TEST_H
#define COSINE_ANGLE_TEST_H

#include <inter/angle/CosineAngle.h>
#include "AngleTestTemplate.h"

class CosineAngleTest : public AngleTestTemplate<CosineAngle>
{
public:

   CosineAngleTest() :
      AngleTestTemplate<CosineAngle>("in/CosineAngle")
   {}

   void testSetUp()
   {  printMethod(TEST_FUNC); } 

   void testForce()
   {
      printMethod(TEST_FUNC);
      b1_ = Vector( 0.2, 0.90,  0.3);
      b2_ = Vector(-0.1, 0.85, -0.4);
      forceTest();
   }

};

TEST_BEGIN(CosineAngleTest)
TEST_ADD(CosineAngleTest, testSetUp)
TEST_ADD(CosineAngleTest, testForce)
TEST_END(CosineAngleTest)

#endif
