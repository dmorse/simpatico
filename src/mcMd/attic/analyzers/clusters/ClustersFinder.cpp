#ifndef MCMD_CLUSTERS_FINDER_CPP
#define MCMD_CLUSTERS_FINDER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ClustersFinder.h"
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <simp/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/boundary/Boundary.h>
#include <util/misc/FileMaster.h>        
#include <util/archives/Serializable_includes.h>

#include <util/format/Dbl.h>
#include <util/misc/ioUtil.h>
#include <sstream>

namespace McMd
{

   using namespace Util;

   /// Constructor.
   ClustersFinder::ClustersFinder(System& system) 
    : SystemAnalyzer<System>(system),
      outputFile_(),
      speciesId_(),
      speciesPtr_(),
      coreId_(),
      cellList_(),
      clusters_(),
      clusterLengths_(),
      histMin_(0),
      histMax_(),
      hist_(),
      nSample_(0),
      isInitialized_(false)
   {  setClassName("ClustersFinder"); }

   /*
   * Clear accumulators.
   */
   void ClustersFinder::setup() 
   {  int nMolecule = system().nMolecule(speciesId_);
      int nAtom = nMolecule * speciesPtr_->nAtom();
      cellList_.allocate(nAtom, system().boundary(), cutoff_);
      nSample_ = 0;
      clusters_.allocate(nMolecule);
      isInitialized_= true;
      if (!isInitialized_) UTIL_THROW("Object is not initialized");
   }

   /*
   * Clear accumulators.
   */
   void ClustersFinder::findClusters(Molecule* molPtr, int clusterId) 
   {
      GArray<int> mNeighbors;
      mNeighbors.reserve(25);
      
      CellList::NeighborArray aNeighbors;

      Molecule::AtomIterator atomIter;
      for (molPtr->begin(atomIter); atomIter.notEnd(); ++atomIter) {                   // Identifying the neighbours! 
         
          if (atomIter->typeId() == coreId_) {                                         // Checks the atomType to make sure it has the right Type.
            cellList_.getNeighbors(atomIter->position(), aNeighbors);                  // Takes neighbors of molecules atom out of cellList.
               //std::cout<<aNeighbors.size()<<"\n";   
        for (int i = 0; i < aNeighbors.size(); i++) {
                   //std::cout<<clusters_[aNeighbors[i]->molecule().id()].clusterId_<<"\n";
                if (clusters_[aNeighbors[i]->molecule().id()].clusterId_ == -1) {
                   //std::cout<<aNeighbors[i]->molecule().id()<<"\t"<<clusterId<<"\n";
                   clusters_[aNeighbors[i]->molecule().id()].clusterId_ = clusterId; 
                   //std::cout<<aNeighbors[i]->molecule().id()<<"\t"<<clusterId<<"\n";
                   mNeighbors.append(aNeighbors[i]->molecule().id());               // Adds neighbor molecule to list of neighbors.
                } else if (clusters_[aNeighbors[i]->molecule().id()].clusterId_ != clusterId) UTIL_THROW("Cluster Clash!"); 
            } 
          }
      }
      
      for (int i = 0; i < mNeighbors.size(); ++i) {
          findClusters(clusters_[mNeighbors[i]].self_, clusterId);
      }
      
      aNeighbors.clear();
      mNeighbors.deallocate();
   }

   /* 
   * Roots out all the Clusters.
   */
   int ClustersFinder::clusterId(Molecule* molPtr)
   {  std::cout << "test I don't even know";
      return  clusters_[molPtr->id()].clusterId_;
   }

   /* 
   * Evaluate end-to-end vectors of all chains, add to ensemble.
   */
   void ClustersFinder::sample(long iStep) 
   {  std::cout << "sample error";
      if (isAtInterval(iStep)) {
         
         int nMolecule = system().nMolecule(speciesId_);
         int clusterId = 0;
         cellList_.clear();
         clusterLengths_.clear();
         hist_.clear();

         for (int i = 0; i < nMolecule; ++i) { 
             clusters_[i].self_ = 0;
             clusters_[i].clusterId_ = -1;
         }

         System::MoleculeIterator molIter;                                             // Loading cellList with atoms.  
         Molecule::AtomIterator atomIter;           
         for (system().begin(speciesId_, molIter); molIter.notEnd(); ++molIter) {
             clusters_[molIter->id()].self_ = molIter.get();
             for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
                 if (atomIter->typeId() == coreId_) {
                    //std::cout<<"\n";
                    system().boundary().shift(atomIter->position());
                    cellList_.addAtom(*atomIter);
                 }
             }              // Atom loop.
         }                  // Molecule loop.
         
         for (system().begin(speciesId_, molIter); molIter.notEnd(); ++molIter) {
              //std::cout<<clusters_[molIter->id()].clusterId_<<"\n";
              if (clusters_[molIter->id()].clusterId_ == -1) {
                findClusters(molIter.get(), clusterId);
                clusterId++;
             }
         }
         clusterLengths_.resize(clusterId);
         for (int i = 0; i < clusterId; ++i) {
             clusterLengths_[i] = 0;
         }
         
         for (int i = 0; i < nMolecule; i++) {
             //if ( clusters_[i].clusterId_ == -1 ) UTIL_THROW("Clusterization not completed!");
             ++clusterLengths_[clusters_[i].clusterId_];
         }

      std::string fileName;
      fileName  = outputFileName();
      fileName += toString(iStep);
      fileName += ".clusters";
      fileMaster().openOutputFile(fileName, outputFile_);
      for (int i = 0; i < clusterLengths_.size(); i++) {
          outputFile_<<"Cluster "<< i+1<<"includes molecules:"<< "\n";
          for (system().begin(speciesId_, molIter); molIter.notEnd(); ++molIter) {
              if(clusters_[molIter->id()].clusterId_ == i) outputFile_<< molIter->id()<<"\t";
          }
          outputFile_<<"\n\n\n";
      }
      outputFile_.close();

      fileName  = outputFileName();
      fileName += toString(iStep);
      fileName += ".dat";
      fileMaster().openOutputFile(fileName, outputFile_);
      outputFile_ << "Cluster Id" <<"\t"<< "Number of Molecules" <<"\n";
      
      for (int i = 0; i < clusterLengths_.size(); i++) {
          hist_.sample(clusterLengths_[i]);
          outputFile_<< i+1 << "\t"<< clusterLengths_[i]<< "\n";
      }

      outputFile_.close();

      fileName  = outputFileName();
      fileName += toString(iStep);
      fileName += ".hist";
      fileMaster().openOutputFile(fileName, outputFile_);
      hist_.output(outputFile_);
      outputFile_.close();
      ++nSample_;
      }
   }

   /*
   * Output results to file after simulation is completed.
   */
   void ClustersFinder::output() 
   {  std::cout << "why?";
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_);
      outputFile_.close();
   }

}
#endif 
