#ifndef SIMP_MOCK_AVERAGE_H
#define SIMP_MOCK_AVERAGE_H

#include <simp/analysis/AverageMixIn.h>
#include <util/param/ParamComposite.h>
#include <util/random/Random.h>

using namespace Util;

namespace Simp {

   class MockAverage : public ParamComposite, public AverageMixIn
   {
   private:

      Random random_;
      std::string outputFileName_;
      double value_;
      int interval_;
   
   
   public:
   
      MockAverage(FileMaster& fileMaster) :
        AverageMixIn(fileMaster)
      { setClassName("MockAverage"); }

      void readParameters(std::istream& in) 
      {
         read<int>(in, "interval", interval_);
         read<std::string>(in, "outputFileName", outputFileName_);
         readNSamplePerBlock(in, *this);
         initializeAccumulator();
      }  

      void loadParameters(Serializable::IArchive& ar) 
      {  
         loadNSamplePerBlock(ar, *this); 
      } 

      void save(Serializable::OArchive& ar) 
      {  
         saveNSamplePerBlock(ar); 
      }  

      void setup() 
      {
         // Open data file
         std::string fileName = outputFileName_;
         fileName += ".dat";
         openOutputFile(fileName);

         // Set random number generator seed
         long seed = 0;
         random_.setSeed(seed);
      }  

      void sample(long iStep) 
      {
         compute();
         updateAccumulator(iStep, interval_);
      }
  
      void output() 
      {
         if (outputFile().is_open()) {
            outputFile().close();
         }

         std::string fileName = outputFileName_;
         fileName += ".prm";
         openOutputFile(fileName);
         writeParam(outputFile());
         outputFile().close();
  
         outputAccumulator(outputFileName_); 
      }  


      double interval() const  
      {  return interval_; }

   protected:
 
      void compute() 
      {  value_ = 1.0 + 0.01*random_.gaussian(); } 
   
      double value()
      {  return value_; }

   };
 
}
#endif
