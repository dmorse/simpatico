#ifndef COSINE_SQ_ANGLE_TEST_H
#define COSINE_SQ_ANGLE_TEST_H

#include <inter/angle/CosineSqAngle.h>
#include "AngleTestTemplate.h"

class CosineSqAngleTest : public AngleTestTemplate<CosineSqAngle>
{
public:

   CosineSqAngleTest() :
      AngleTestTemplate<CosineSqAngle>("in/CosineSqAngle")
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

TEST_BEGIN(CosineSqAngleTest)
TEST_ADD(CosineSqAngleTest, testSetUp)
TEST_ADD(CosineSqAngleTest, testForce)
TEST_END(CosineSqAngleTest)

#endif
