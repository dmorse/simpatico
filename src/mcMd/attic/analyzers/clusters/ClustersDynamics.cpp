#ifndef MCMD_CLUSTERS_DYNAMICS_CPP
#define MCMD_CLUSTERS_DYNAMICS_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ClustersDynamics.h"
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
   ClustersDynamics::ClustersDynamics(System& system) 
    : SystemAnalyzer<System>(system),
      outputFile_(),
      oldClusters_(system),
      newClusters_(system),
      oldClustersPtr_(&oldClusters_),
      newClustersPtr_(&newClusters_),
      speciesId_(),
      speciesPtr_(),
      coreId_(),
      criterion_(),
      histMin_(0),
      histMax_(),
      isInitialized_(false)
   {  setClassName("ClustersDynamics"); }

   /// Read parameters from file, and allocate arrays.
   void ClustersDynamics::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);

      read<int>(in, "speciesId", speciesId_);
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      if (speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("speciesId > nSpecies");
      }

      oldClusters_.speciesPtr_ = &system().simulation().species(speciesId_);
      newClusters_.speciesPtr_ = &system().simulation().species(speciesId_);

      oldClusters_.speciesId_ = speciesId_;
      newClusters_.speciesId_ = speciesId_;



      read<int>(in, "coreId", coreId_);
      if (coreId_ < 0) {
         UTIL_THROW("Negative coreId");
      }
      read<double>(in, "criterion", criterion_);
      if (criterion_ < 0) {
         UTIL_THROW("Negative criterion");
      }
      oldClusters_.coreId_ = coreId_;
      newClusters_.coreId_ = coreId_;
      
      read<double>(in, "cutoff", cutoff_);
      if (cutoff_ < 0) {
         UTIL_THROW("Negative cutoff");
      }
      oldClusters_.cutoff_ = cutoff_;
      newClusters_.cutoff_ = cutoff_;

      read<int>(in,"histMin", histMin_);
      read<int>(in,"histMax", histMax_);
      oldClusters_.hist_.setParam(histMin_, histMax_);
      newClusters_.hist_.setParam(histMin_, histMax_);

      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void ClustersDynamics::loadParameters(Serializable::IArchive& ar)
   { 
      // Load (everything but accumulators_)
      Analyzer::loadParameters(ar);

      loadParameter<int>(ar,"speciesId", speciesId_);
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      if (speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("speciesId > nSpecies");
      }

      oldClusters_.speciesPtr_ = &system().simulation().species(speciesId_);
      newClusters_.speciesPtr_ = &system().simulation().species(speciesId_);
      oldClusters_.speciesId_ = speciesId_;
      newClusters_.speciesId_ = speciesId_;

      loadParameter<int>(ar, "coreId", coreId_);
      if (coreId_ < 0) {
         UTIL_THROW("Negative coreId");
      }
      oldClusters_.coreId_ = coreId_;
      newClusters_.coreId_ = coreId_;

      loadParameter<double>(ar, "criterion", criterion_);
      if (criterion_ < 0) {
         UTIL_THROW("Negative criterion");
      }
      loadParameter<double>(ar, "cutoff", cutoff_);
      if (cutoff_ < 0) {
         UTIL_THROW("Negative cutoff");
      }
      oldClusters_.cutoff_ = cutoff_;
      newClusters_.cutoff_ = cutoff_;

      loadParameter<int>(ar,"histMin", histMin_);
      loadParameter<int>(ar,"histMax", histMax_);
      oldClusters_.hist_.setParam(histMin_, histMax_);
      newClusters_.hist_.setParam(histMin_, histMax_);

      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   void ClustersDynamics::save(Serializable::OArchive& ar)
   {  
      ar & *this; }

   /*
   * Clear accumulators.
   */
   void ClustersDynamics::setup() 
   {  oldClusters_.setup();
      newClusters_.setup();
      if (!isInitialized_) UTIL_THROW("Object is not initialized");
   }

   /* 
   * Evaluate end-to-end vectors of all chains, add to ensemble.
   */
   void ClustersDynamics::sample(long iStep) 
   { if (isAtInterval(iStep)) {
         newClustersPtr_->sample(iStep);

         if (iStep == 0) {
            oldClustersPtr_->sample(iStep);
         }         
         int clustersCorrelations_[newClustersPtr_->clusterLengths_.size()][oldClustersPtr_->clusterLengths_.size()];
          
         for (int i = 0; i < newClustersPtr_->clusterLengths_.size(); i++) {
             for (int j = 0; j < oldClustersPtr_->clusterLengths_.size(); j++) {
                 clustersCorrelations_[i][j] = 0;
             }
         }

         System::MoleculeIterator molIter;
         for (system().begin(speciesId_, molIter); molIter.notEnd(); ++molIter) {
             clustersCorrelations_[newClustersPtr_->clusters_[molIter->id()].clusterId_][oldClustersPtr_->clusters_[molIter->id()].clusterId_] += 1 ;
         }

         bool usedIds_[newClustersPtr_->clusterLengths_.size()];
         int idConversions_[newClustersPtr_->clusterLengths_.size()];

         for (int i = 0; i < newClustersPtr_->clusterLengths_.size(); i++) {
             usedIds_[i] = false;
             idConversions_[i] = 0;
         }
         int id_ = 0;
         for (int i = 0; i < newClustersPtr_->clusterLengths_.size(); i++) {
           for (int j = 0; j < oldClustersPtr_->clusterLengths_.size(); j++) {
             //this only works if the criterion > 50%
                 if (clustersCorrelations_[i][j] > criterion_ * oldClustersPtr_->clusterLengths_[i]) id_ = i;
             }
             usedIds_[id_] = true;
             idConversions_[i] = id_;
         }
         int usedIdCounter_ = 0;
         for (int i = 0; i < newClustersPtr_->clusterLengths_.size(); i++) {
             if (idConversions_[i] == 0) {
                while (usedIds_[usedIdCounter_] == true) usedIdCounter_++;
                idConversions_[i] = usedIdCounter_;
                usedIds_[usedIdCounter_] = true;
             }  
         }

/*         for (system().begin(speciesId_, molIter); molIter.notEnd(); ++molIter) {
             newClustersPtr_->clusters_[molIter->id()].clusterId_ = idConversions_[newClustersPtr_->clusters_[molIter->id()].clusterId_];
         } */

         std::string fileName;
         fileName  = outputFileName();
         fileName += toString(iStep);
         fileName += ".dynamics";
         fileMaster().openOutputFile(fileName, outputFile_);
         for (int i = 0; i < newClustersPtr_->clusterLengths_.size(); i++) {
             outputFile_<<"<"<< i <<">"<<"("<<newClustersPtr_->clusterLengths_[i]<<") = ";
             for (int j=0; j < oldClustersPtr_->clusterLengths_.size(); j++) {
                  if (clustersCorrelations_[i][j] != 0){
                    outputFile_<<"<"<< j <<">["<<oldClustersPtr_->clusterLengths_[j]<<"]("<<clustersCorrelations_[i][j]<<") ";
                 }
             }
             outputFile_<<"\n";
         }
         outputFile_.close();

         ClustersFinder *tmp = oldClustersPtr_;
         oldClustersPtr_ = &(*newClustersPtr_);
         newClustersPtr_ = &(*tmp);
      }
   }

   /*
   * Output results to file after simulation is completed.
   */
   void ClustersDynamics::output() 
   {  std::cout << "testing6";
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_);
      outputFile_.close();
   }

}
#endif 
