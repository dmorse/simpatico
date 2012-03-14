#ifndef COSINE_DIHEDRAL_TEST_H
#define COSINE_DIHEDRAL_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <inter/dihedral/CosineDihedral.h>
#include <mcMd/chemistry/Atom.h>
#include <util/random/Random.h>
#include <util/containers/RArray.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Inter;

class CosineDihedralTest : public UnitTest 
{

private:

   CosineDihedral dihedralPotential;
   Vector        r[3];

public:

   void setUp()
   {
      // Set Boundary Lengths
      dihedralPotential.setNDihedralType(1);

      // Read parameters from file
      std::ifstream in("dihedral/in/CosineDihedral");
      dihedralPotential.readParam(in);
      in >> r[0];
      in >> r[1];
      in >> r[2];
      in.close();

      //r[0] = Vector(1.0, 0.3, 0.5);
      //r[1] = Vector(1.0, 1.0, 1.0);
      //r[2] = Vector(1.0, 1.0, 0.4);
   }


   void tearDown()
   {}


   void testSetUp() 
   {
      printMethod(TEST_FUNC);
   }


   void testWrite() {
      printMethod(TEST_FUNC);

      // Verbose output
      if (verbose() > 0) { 
         dihedralPotential.writeParam(std::cout);
      }
   }


   void testEnergy() 
   {
      Vector re[3], force;
      double energy, eps, newEnergy;
     
      printMethod(TEST_FUNC);
      std::cout << std::endl;
      
      energy = dihedralPotential.energy(r[0], r[1], r[2], 0);
      std::cout << r[0] << std::endl;
      std::cout << r[1] << std::endl;
      std::cout << r[2] << std::endl;
      std::cout << "dihedralPotential.energy = " << energy << std::endl;

      eps = 0.00001;
      re[0] = Vector(eps, 0.0, 0.0);
      re[1] = Vector(0.0, eps, 0.0);
      re[2] = Vector(0.0, 0.0, eps);

      std::cout << "dihedralPotential.energy perturbation: " << std::endl;
      for (int j = 0; j < 3; ++j) { // bond vectors
         for (int k = 0; k < 3; ++k) { // components
            r[j] += re[k];
            newEnergy = dihedralPotential.energy(r[0], r[1], r[2], 0);
            r[j] -= re[k];
            force[k] = (newEnergy - energy) / eps;
         }
         std::cout << force << std::endl;
      }
   }

   void testForce() 
   {
      Vector f1, f2, f3;
 
      printMethod(TEST_FUNC);
      std::cout << std::endl;
     
      dihedralPotential.force(r[0], r[1], r[2], f1, f2, f3, 0);

      std::cout << "dihedralPotential.force: " << std::endl;
      std::cout << f1 << std::endl;
      std::cout << f2 << std::endl;
      std::cout << f3 << std::endl;
   }

};

TEST_BEGIN(CosineDihedralTest)
TEST_ADD(CosineDihedralTest, testSetUp)
TEST_ADD(CosineDihedralTest, testWrite)
TEST_ADD(CosineDihedralTest, testEnergy)
TEST_ADD(CosineDihedralTest, testForce)
TEST_END(CosineDihedralTest)

#endif
