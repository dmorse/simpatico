#ifndef SIMP_MOCK_AVERAGE_LIST_H
#define SIMP_MOCK_AVERAGE_LIST_H

#include <simp/analysis/AverageListMixIn.h>
#include <util/param/ParamComposite.h>
#include <util/random/Random.h>

using namespace Util;

namespace Simp {

   class MockAverageList : public ParamComposite, public AverageListMixIn
   {
   private:

      Random random_;
      std::string outputFileName_;
      int interval_;
      int nValue_;
   
   public:
   
      MockAverageList(FileMaster& fileMaster) :
        AverageListMixIn(fileMaster),
        outputFileName_("data"),
        interval_(1),
        nValue_(2)
      { setClassName("MockAverageList"); }

      void readParameters(std::istream& in) 
      {
         read<int>(in, "interval", interval_);
         read<std::string>(in, "outputFileName", outputFileName_);
         readNSamplePerBlock(in, *this);
         initializeAccumulators(nValue_);
         setName(0, "A");
         setName(1, "B");
      }  

      void loadParameters(Serializable::IArchive& ar) 
      {  
         loadParameter(ar, "interval", interval_);
         loadParameter(ar, "outputFileName", outputFileName_);
         loadNSamplePerBlock(ar, *this); 
         loadAccumulators(ar);
      } 

      void save(Serializable::OArchive& ar) 
      {  
         ar << interval_;
         ar << outputFileName_;
         saveNSamplePerBlock(ar); 
         saveAccumulators(ar); 
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
         updateAccumulators(iStep, interval_);
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
  
         outputAccumulators(outputFileName_); 
      }  

      double interval() const  
      {  return interval_; }

   protected:
 
      void compute() 
      {  
         double average, current;
         for (int i = 0; i < nValue_; ++i) {
            average = (double)i;
            current =  average + 0.01*random_.gaussian(); 
            setValue(i, current);
         } 
      }
   
   };
 
}
#endif
