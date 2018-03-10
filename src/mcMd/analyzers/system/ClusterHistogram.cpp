#ifndef MCMD_CLUSTER_HISTOGRAM_CPP
#define MCMD_CLUSTER_HISTOGRAM_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ClusterHistogram.h"
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>

#include <simp/boundary/Boundary.h>
#include <simp/species/Species.h>

#include <util/misc/FileMaster.h>        
#include <util/archives/Serializable_includes.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/misc/ioUtil.h>
#include <util/space/Tensor.h>
#include <util/containers/DArray.h>
#include <sstream>

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /// Constructor.
   ClusterHistogram::ClusterHistogram(System& system) 
    : SystemAnalyzer<System>(system),
      identifier_(system),
      hist_(),
      outputFile_(),
      speciesId_(),
      atomTypeId_(),
      cutoff_(),
      histMin_(),
      histMax_(),
      nSample_(0),
      isInitialized_(false)
   {  setClassName("ClusterHistogram"); }

   /// Read parameters from file, and allocate arrays.
   void ClusterHistogram::readParameters(std::istream& in) 
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

      read<int>(in, "atomTypeId", atomTypeId_);
      if (atomTypeId_ < 0) {
         UTIL_THROW("Negative atomTypeId");
      }
      if (atomTypeId_ >= system().simulation().nAtomType()) {
         UTIL_THROW("nTypeId >= nAtomType");
      }

      read<double>(in, "cutoff", cutoff_);
      if (cutoff_ < 0) {
         UTIL_THROW("Negative cutoff");
      }

      // Initialize ClusterIdentifier
      identifier_.initialize(speciesId_, atomTypeId_, cutoff_);
      read<int>(in,"histMin", histMin_);
      read<int>(in,"histMax", histMax_);
      hist_.setParam(histMin_, histMax_);
      hist_.clear();
    
      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void ClusterHistogram::loadParameters(Serializable::IArchive& ar)
   {
      // Load interval and outputFileName
      Analyzer::loadParameters(ar);

      loadParameter<int>(ar,"speciesId", speciesId_);
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      if (speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("speciesId > nSpecies");
      }

      loadParameter<int>(ar, "atomTypeId", atomTypeId_);
      if (atomTypeId_ < 0) {
         UTIL_THROW("Negative atomTypeId");
      }

      loadParameter<double>(ar, "cutoff", cutoff_);
      if (cutoff_ < 0) {
         UTIL_THROW("Negative cutoff");
      }

      identifier_.initialize(speciesId_, atomTypeId_, cutoff_);
      loadParameter<int>(ar, "histMin", histMin_);
      loadParameter<int>(ar, "histMax", histMax_);
      ar >> hist_;

      ar >> nSample_;

      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   void ClusterHistogram::save(Serializable::OArchive& ar)
   {
      Analyzer::save(ar);
      ar & speciesId_;
      ar & atomTypeId_;
      ar & cutoff_;
      ar & histMin_;
      ar & histMax_;
      ar & hist_;
      ar & nSample_;
   }

   /*
   * Clear accumulators.
   */
   void ClusterHistogram::setup() 
   {  
      if (!isInitialized_) UTIL_THROW("Object is not initialized");
      hist_.clear();
      nSample_ = 0;
   }

   /* 
   * Sample data by calling ClusterIdentifier::identifyClusters.
   */
   void ClusterHistogram::sample(long iStep) 
   { 
      if (isAtInterval(iStep)) {
         //Identifies all clusters
         identifier_.identifyClusters();
         //Adds each cluster to the histogram of all clusters
         for (int i = 0; i < identifier_.nCluster(); i++) {
             hist_.sample(identifier_.cluster(i).size());
         }
         ++nSample_;
         fileMaster().openOutputFile(outputFileName(".clusters"+toString(iStep)),outputFile_);
         //Writes all of the clusters and their component molecules
         Cluster thisCluster;
         ClusterLink* thisClusterStart;
         ClusterLink* next;
         Molecule thisMolecule;
         //Loop over each cluster
         for (int i = 0; i < identifier_.nCluster(); i++) {
             thisCluster = identifier_.cluster(i);
             thisClusterStart = thisCluster.head();
             outputFile_ << i << "	" ;
             //List out every molecule in that cluster
             while (thisClusterStart) {
                next = thisClusterStart->next();
                thisMolecule = thisClusterStart->molecule();
                outputFile_ << thisMolecule.id() << "  ";
                thisClusterStart = next;
             }
             outputFile_ << "\n";
         }
         outputFile_.close();

         //Calculate, store, and write the micelle center of mass
         fileMaster().openOutputFile(outputFileName(".COMs"+toString(iStep)),outputFile_);
         //comArray;
         Vector clusterCOM;
         Vector r0;
         Vector dr;
         Tensor moment;
         Tensor rgDyad;
         DArray<Vector> allCOMs;
         DArray<Tensor> allMoments;
         allCOMs.allocate(identifier_.nCluster());
         allMoments.allocate(identifier_.nCluster());
         Molecule::ConstAtomIterator atomIter;
         for (int i = 0; i < identifier_.nCluster(); i++) {
             thisCluster = identifier_.cluster(i);
             outputFile_ << i << "	" ;
             //For that cluster, calculate the center of mass
             clusterCOM = thisCluster.clusterCOM(atomTypeId_, system().boundary());
             outputFile_ << clusterCOM;
             outputFile_ << "\n";
             allCOMs[i] = clusterCOM;
             //Calculate Rg
             moment = thisCluster.momentTensor(atomTypeId_, system().boundary());
             allMoments[i] = moment;
         }
         outputFile_.close();
         fileMaster().openOutputFile(outputFileName(".momentTensors"+toString(iStep)),outputFile_);
         for (int i = 0; i < identifier_.nCluster(); i++) {
             outputFile_ << i << "	" << allMoments[i] << "\n";
           
         }
         outputFile_.close();
      }
   }

   /*
   * Output results to file after simulation is completed.
   */
   void ClusterHistogram::output() 
   {
      // Write parameter file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_);
      outputFile_.close();

      // Write histogram output
      fileMaster().openOutputFile(outputFileName(".hist"), outputFile_);
      // hist_.output(outputFile_);
      int min = hist_.min();
      int nBin = hist_.nBin();
      for (int i = 0; i < nBin; ++i) {
         outputFile_ << Int(i + min) << "  " 
                     <<  Dbl(double(hist_.data()[i])/double(nSample_)) << "\n";
      }
      outputFile_.close();
   }

}
#endif 
