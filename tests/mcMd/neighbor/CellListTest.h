#ifndef MCMD_CELL_LIST_TEST_H
#define MCMD_CELL_LIST_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <mcMd/neighbor/CellList.h>
#include <util/boundary/Boundary.h>
#include <mcMd/chemistry/Atom.h>
#include <util/space/Vector.h>
#include <util/random/Random.h>
#include <util/containers/RArray.h>

#include <iostream>

using namespace Util;
using namespace McMd;


class CellListTest : public UnitTest 
{

private:

   CellList cellList;
   Boundary boundary;

public:

   void setUp()
   {}

   void tearDown()
   {}

   void testMakeGrid()
   {
      printMethod(TEST_FUNC);
      double cutoff  = 1.2;

      Vector  Lin;
      Lin[0] =  2.0;
      Lin[1] =  3.0;
      Lin[2] =  4.0;
      boundary.setOrthorhombic(Lin);  
      cellList.makeGrid(boundary, cutoff);

      TEST_ASSERT(cellList.gridDimension(0) == 1);
      TEST_ASSERT(cellList.gridDimension(1) == 2);
      TEST_ASSERT(cellList.gridDimension(2) == 3);

      //printf("\n");
      //printf("CellList.numCells: %i %i %i \n", 
      //cellList.gridDimension(0),cellList.gridDimension(1),cellList.gridDimension(2));
   }

   void testShiftCellCoordAxis(){
      int x, y, z;
      double cutoff  = 1.2;
      printMethod(TEST_FUNC);
      //printf("CellListTest::testShiftCellCoordAxis\n");

      Vector  Lin;
      Lin[0] =  2.0;  // 1 cell
      Lin[1] =  3.0;  // 2 cells
      Lin[2] =  4.0;  // 3 cells
      boundary.setOrthorhombic(Lin);  
      cellList.allocate(10, boundary, cutoff);
      cellList.makeGrid(boundary, cutoff);

      x =  1;
      y =  3;
      z = -2;
      TEST_ASSERT(cellList.shiftCellCoordAxis(0,x) == 0);
      TEST_ASSERT(cellList.shiftCellCoordAxis(1,y) == 1);
      TEST_ASSERT(cellList.shiftCellCoordAxis(2,z) == 1);

   }

   void testCellIndexCoord(){
      int i, x, y, z, xn, yn, zn;
      double cutoff  = 1.2;
      printMethod(TEST_FUNC);

      Vector Lin;
      Lin[0] =  2.0;
      Lin[1] =  3.0;
      Lin[2] =  4.0;
      boundary.setOrthorhombic(Lin);  
      cellList.allocate(10, boundary, cutoff);
      cellList.makeGrid(boundary, cutoff);

      // Test consistency of CellIndexFromCoord and CellCoordFromIndex
      x =  0;
      y =  1;
      z =  2;
      i = cellList.cellIndexFromCoord(x,y,z);
      cellList.cellCoordFromIndex(i, xn,yn,zn);

      TEST_ASSERT(xn == x);
      TEST_ASSERT(yn == y);
      TEST_ASSERT(zn == z);
      TEST_ASSERT(i == 5);

      // printf("\n");
      // printf("x, y, z, i: %i %i %i %i \n", x, y, z, i);

   }

   void testInitialize(){
      double cutoff  = 1.2;
      printMethod(TEST_FUNC);

      Vector Lin;
      Lin[0] = 2.0;
      Lin[1] = 3.0;
      Lin[2] = 4.0;
      boundary.setOrthorhombic(Lin);  
      cellList.allocate(10, boundary, cutoff);
      cellList.makeGrid(boundary, cutoff);

      TEST_ASSERT(cellList.totCells_ == 6);
      TEST_ASSERT(cellList.gridDimension(0) == 1);
      TEST_ASSERT(cellList.gridDimension(1) == 2);
      TEST_ASSERT(cellList.gridDimension(2) == 3);
      TEST_ASSERT(eq(cellList.invCellWidths_[0], 1.0));
      TEST_ASSERT(eq(cellList.invCellWidths_[1], 2.0));
      TEST_ASSERT(eq(cellList.invCellWidths_[2], 3.0));

   }

   void testCellIndexFromPosition(){
      Vector r;
      double cutoff  = 1.2;
      int    i;
      printMethod(TEST_FUNC);

      Vector Lin;
      Lin[0] =  2.8; // 2 cell
      Lin[1] =  3.0; // 2 cells
      Lin[2] =  4.0; // 3 cells
      boundary.setOrthorhombic(Lin);  
      cellList.allocate(10, boundary, cutoff);
      cellList.makeGrid(boundary, cutoff);

      TEST_ASSERT(cellList.gridDimension(0) == 2);
      TEST_ASSERT(cellList.gridDimension(1) == 2);
      TEST_ASSERT(cellList.gridDimension(2) == 3);
      TEST_ASSERT(cellList.YZCells_ == 6);

      r = Vector(1.3, 1.3, 2.5);
       // x=0, y=0, z=1  
      i = cellList.cellIndexFromPosition(r);
      // Assert: i = 0*6 + 0*3 + 1 = 1
      TEST_ASSERT(i == 1); 

      r = Vector(1.45, 1.3, 3.5);
      // x=1, y=0, z=2
      i = cellList.cellIndexFromPosition(r);
      // Assert: i= 1*6 + 0*3 + 2 = 8
      TEST_ASSERT(i == 8);   

   }


   void testAddAtom(){
      const int    maxNAtom = 5;
      const double cutoff  = 1.2;

      printMethod(TEST_FUNC);

      // Create a Boundary
      Vector Lin(2.0, 3.3, 2.5);
      //Lin[0] =  2.0;  // 1 cell,  x=0
      //Lin[1] =  3.3;  // 2 cells, y=0
      //Lin[2] =  2.5;  // 2 cells, z=1, icell=1
      boundary.setOrthorhombic(Lin);  
      cellList.allocate(maxNAtom, boundary, cutoff);
      cellList.makeGrid(boundary, cutoff);

      // Allocate and initialize an array of Atoms
      RArray<Atom> atoms;
      Atom::allocate(maxNAtom, atoms);

      Vector position(1.0, 1.3, 1.3);
      int    atomId = 3;
      atoms[atomId].position() = position;
      int cellId = cellList.cellIndexFromPosition(atoms[atomId].position());
      TEST_ASSERT(1 == cellId);

      cellList.addAtom(atoms[atomId]);

      TEST_ASSERT(cellList.cells_[cellId].nAtomCell_ == 1);
      TEST_ASSERT(cellList.cells_[cellId].firstEmptyPos_ == 1);

      TEST_ASSERT(cellList.cellTags_[atomId].cellId == cellId);
      TEST_ASSERT(cellList.cellTags_[atomId].cellPos == 0);

      try { 
         cellList.isValid(1);
      } catch (Exception e) { 
          e.write(std::cout);
          TEST_ASSERT(0);
      }

      Atom::deallocate();
   }

   void testDeleteAtom()
   {
      printMethod(TEST_FUNC);

      const int    maxNAtom = 10;
      const double cutoff   = 1.2;

      // Create a OrthorhombicBoundary
      Vector Lin(2.0, 3.0, 2.5);
      //Lin[0] =  2.0;  // 1 cell,  x = 0
      //Lin[1] =  3.0;  // 2 cells, y = 1
      //Lin[2] =  2.5;  // 2 cells, z = 1
      boundary.setOrthorhombic(Lin);  
      cellList.allocate(maxNAtom, boundary, cutoff);
      cellList.makeGrid(boundary, cutoff);

      // Allocate and initialize an array of Atoms
      RArray<Atom>  atoms;
      Atom::allocate(maxNAtom, atoms);

      Vector r(1.0, 1.6, 1.3);
      int cellId = cellList.cellIndexFromPosition(r);
      TEST_ASSERT(cellId == 3);

      int    atomId = 3;
      atoms[atomId].position() = r;
      cellList.addAtom(atoms[atomId]);
      cellList.deleteAtom(atoms[atomId]);

      TEST_ASSERT(cellList.cells_[cellId].nAtomCell_ == 0);
      TEST_ASSERT(cellList.cells_[cellId].firstEmptyPos_ == 0);

      try { 
         cellList.isValid(0);
      } catch (Exception e) { 
          e.write(std::cout);
          TEST_ASSERT(0);
      }

      Atom::deallocate();
   }


   void testAddDeleteAtoms(){
      // Add N_PART atoms to the same cell
      double cutoff  = 1.2;

      #define N_PART   2
      Vector  r[N_PART] = { 
                             Vector(1.0, 1.6, 1.0), 
                             Vector(0.5, 1.8, 1.2) 
                          };
      int cellId[N_PART];
      int i,k;
      int j=0; // Index of atom to be removed

      printMethod(TEST_FUNC);

      // Initialize an OrthorhombicBoundary and CellList geometry
      Vector Lin(2.0, 3.0, 4.0);
      //Lin[0] =  2.0;  // 1 cell,   x1=0 x2=0
      //Lin[1] =  3.0;  // 2 cells,  x1=2 x2=2
      //Lin[2] =  4.0;  // 3 cells,  x1=0 x2=0
      boundary.setOrthorhombic(Lin);  
      cellList.allocate(N_PART, boundary, cutoff);
      cellList.makeGrid(boundary, cutoff);

      // Allocate and initialize an array of N_PART atoms
      RArray<Atom>  atoms;
      Atom::allocate(N_PART, atoms);

      // Add N_PART atoms to the same cell
      for (i=0; i < N_PART; i++){
         cellId[i] = cellList.cellIndexFromPosition(r[i]);
         atoms[i].position() = r[i];
         cellList.addAtom(atoms[i]);

         TEST_ASSERT(cellList.cells_[cellId[i]].nAtomCell_ == i+1);
         TEST_ASSERT(cellList.cells_[cellId[i]].firstEmptyPos_ == i+1);
         TEST_ASSERT(cellList.cells_[cellId[i]].atoms_[i] == &atoms[i]);

         TEST_ASSERT(cellList.cellTags_[i].cellId  == cellId[i]);
         TEST_ASSERT(cellList.cellTags_[i].cellPos == i);
      }

      i = cellId[0];

      TEST_ASSERT(cellList.cells_[i].nAtomCell_ == N_PART);
      TEST_ASSERT(cellList.cells_[i].firstEmptyPos_ == N_PART);

      // Remove atom j from cell (not the last one)
      cellList.deleteAtom(atoms[j]);

      TEST_ASSERT(cellList.cells_[i].nAtomCell_ == N_PART - 1);
      TEST_ASSERT(cellList.cells_[i].firstEmptyPos_ == j);

      // Put back atom j, with modified position and label
      r[j] = Vector(0.9, 2.4, 1.1);
      atoms[j].position() = r[j];
      i = cellList.cellIndexFromPosition(atoms[j].position());
      TEST_ASSERT(i == cellId[j]);  // Check that its the same cell
      cellList.addAtom(atoms[j]);

      TEST_ASSERT(cellList.cells_[i].nAtomCell_ == N_PART);
      TEST_ASSERT(cellList.cells_[i].firstEmptyPos_ == N_PART);

      for (k=0; k < N_PART; k++) {

         TEST_ASSERT(cellList.cells_[i].atoms_[k] == &atoms[k]);

         TEST_ASSERT(cellList.cellTags_[k].cellId  == i);
         TEST_ASSERT(cellList.cellTags_[k].cellPos == k);
      }

      for (k=N_PART; k < Cell::MaxAtomCell; k++) {
         TEST_ASSERT(cellList.cells_[i].atoms_[k] == 0);
      }

      TEST_ASSERT(cellList.cellTags_[j].cellId  == i);
      TEST_ASSERT(cellList.cellTags_[j].cellPos == j);

      try { 
         cellList.isValid(N_PART);
      } catch (Exception e) { 
          e.write(std::cout);
          TEST_ASSERT(0);
      }

      #undef N_PART

      Atom::deallocate();
   }


   void testAddRandomAtoms(){
      const int nAtom = 20;
      Vector    pos;
      Random    random;
      double    cutoff  = 1.2;
      int       i;

      printMethod(TEST_FUNC);

      // Setup CellList
      Vector Lin(2.0, 3.0, 4.0);
      boundary.setOrthorhombic(Lin);  
      cellList.allocate(nAtom, boundary, cutoff);
      cellList.makeGrid(boundary, cutoff);

      // Allocate Atom 
      RArray<Atom>  atoms;
      Atom::allocate(nAtom, atoms);

      // Create new random number generator
      random.setSeed(1098640);

      // Place atoms at random
      for (i=0; i < nAtom; i++) {
         boundary.randomPosition(random, pos);
         atoms[i].setTypeId(1);
         atoms[i].position() = pos;

         try {
            cellList.addAtom(atoms[i]);
         } 
         catch (Exception e) {
            e.write(std::cout);
            TEST_ASSERT(0);
         }

      }

      try {
         cellList.isValid(nAtom);
      }
      catch (Exception e) {
         e.write(std::cout);
         TEST_ASSERT(0);
      }

      Atom::deallocate();
   }


   void testBuild()
   {
      printMethod(TEST_FUNC);
      const int nAtom = 20;
      Vector       pos;
      Random       random;
      double       cutoff  = 1.2;
      int          i;

      // Setup CellList
      Vector Lin(2.0, 3.0, 4.0);
      boundary.setOrthorhombic(Lin);  
      cellList.allocate(nAtom, boundary, cutoff);
      cellList.makeGrid(boundary, cutoff);

      // Allocate Atom 
      RArray<Atom> atoms;
      Atom::allocate(nAtom, atoms);

      // Create new random number generator
      random.setSeed(1098640);

      // Place atoms at random
      for (i=0; i < nAtom; i++) {
         boundary.randomPosition(random, pos);
         atoms[i].setTypeId(1);
         atoms[i].position() = pos;
      }

      // Add all atoms to the cell list
      cellList.clear();
      for (i=0; i < nAtom; i++) {
         try {
            cellList.addAtom(atoms[i]);
         } 
         catch (Exception e) {
            e.write(std::cout);
            TEST_ASSERT(0);
         }
      }

      // Check validity
      try {
         cellList.isValid(nAtom);
      }
      catch (Exception e) {
         e.write(std::cout);
         TEST_ASSERT(0);
      }

      Atom::deallocate();
   }


   void testUpdateAtomCell(){
      const int nAtom = 20;
      Vector       pos;
      Random       random;
      double       cutoff  = 1.2;
      int          i;
      //double     x, y, z;

      printMethod(TEST_FUNC);

      // Setup CellList
      Vector Lin(2.0, 3.0, 4.0);
      boundary.setOrthorhombic(Lin);  
      cellList.allocate(nAtom, boundary, cutoff);
      cellList.makeGrid(boundary, cutoff);

      // Allocate atoms array
      RArray<Atom> atoms;
      Atom::allocate(nAtom, atoms);

      // Create new random number generator
      random.setSeed(1098640);

      // Place atoms at random
      for (i=0; i < nAtom; ++i) {
         boundary.randomPosition(random, pos);
         atoms[i].setTypeId(1);
         atoms[i].position() = pos;

         try {
            cellList.addAtom(atoms[i]);
         } 
         catch (Exception e) {
            e.write(std::cout);
            TEST_ASSERT(0);
         }
      }

      // Move every atoms to a new random location
      for (i=0; i < nAtom; ++i) {
         boundary.randomPosition(random, pos);
         try {
            cellList.updateAtomCell(atoms[i], pos);
         }
         catch (Exception e) {
            e.write(std::cout);
            TEST_ASSERT(0);
         }
      }


      try {
         cellList.isValid(nAtom);
      }
      catch (Exception e) {
         e.write(std::cout);
         TEST_ASSERT(0);
      }

      Atom::deallocate();
   }


   void testGetNeighbors()
   {
      printMethod(TEST_FUNC);
      const int  nAtom = 20;
      double     cutoff  = 1.2;
      int        i;

      // Set up CellList
      Vector Lin(2.0, 3.0, 4.0);
      boundary.setOrthorhombic(Lin);  
      cellList.allocate(nAtom, boundary, cutoff);
      cellList.makeGrid(boundary, cutoff);

      // Allocate Atoms and initialize Ids
      RArray<Atom>  atoms;
      Atom::allocate(nAtom, atoms);

      // Place Atoms at random
      Vector        pos;
      Random        random;
      random.setSeed(1098640);
      for (i=0; i < nAtom; i++) {
         boundary.randomPosition(random, pos);
         atoms[i].setTypeId(1);
         atoms[i].position() = pos;
      }

      // Add all atoms to the cell list
      cellList.clear();
      for (i=0; i < nAtom; i++) {
         cellList.addAtom(atoms[i]);
      }

      // Find Cell neighbors of a Cell
      CellList::NeighborArray neighborPtrs;
      int           ic, nNeighbors, nInCell;
      //neighborPtrs.allocate(nAtom);
      ic = 5;
      cellList.getCellNeighbors(ic, neighborPtrs, nInCell);
      nNeighbors = neighborPtrs.size();

      std::cout << std::endl;
      std::cout << "Cell Index         = " << ic << std::endl;
      std::cout << "Number in cell      = " << nInCell    << std::endl;
      std::cout << "Number of neighbors = " << nNeighbors << std::endl;

      //const Atom *atomPtr;
      //for (i=0; i < nNeighbors; i++) {
      //   atomPtr = neighborPtrs[i];
      //   std::cout << " Id = " << atomPtr->id() << 
      //           " Position = " << atomPtr->position() << std::endl;
      //}
      
      Vector r(0.9, 2.0, 3.0);
      cellList.getNeighbors(r, neighborPtrs);
      nNeighbors = neighborPtrs.size();
      std::cout << "Cell Index         = " 
           << cellList.cellIndexFromPosition(r) << std::endl;
      std::cout << "Number of neighbors = " << nNeighbors << std::endl;

      //for (i=0; i < nNeighbors; i++) {
      //   atomPtr = neighborPtrs[i];
      //   std::cout << " Id = " << atomPtr->id() << 
      //           " Position = " << atomPtr->position() << std::endl;
      //}

      Atom::deallocate();
   }

   void writeCellConfiguration()
   {
      printf("numCells: %i %i %i \n", 
          cellList.gridDimension(0), cellList.gridDimension(1), cellList.gridDimension(2));
      printf("YZCells, totCells: %i %i \n", 
          cellList.YZCells_, cellList.totCells_);
      printf("cellWidth: %10f %10f %10f \n", 
          cellList.invCellWidths_[0],cellList.invCellWidths_[1],cellList.invCellWidths_[2]);
   }

};

TEST_BEGIN(CellListTest)
TEST_ADD(CellListTest, testMakeGrid)
TEST_ADD(CellListTest, testShiftCellCoordAxis)
TEST_ADD(CellListTest, testCellIndexCoord)
TEST_ADD(CellListTest, testInitialize)
TEST_ADD(CellListTest, testCellIndexFromPosition)
TEST_ADD(CellListTest, testAddAtom)
TEST_ADD(CellListTest, testDeleteAtom)
TEST_ADD(CellListTest, testAddDeleteAtoms)
TEST_ADD(CellListTest, testAddRandomAtoms)
TEST_ADD(CellListTest, testBuild)
TEST_ADD(CellListTest, testUpdateAtomCell)
TEST_ADD(CellListTest, testGetNeighbors)
TEST_END(CellListTest)

#endif
