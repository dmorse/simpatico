#include "LabelBuffer.h"
#include <util/global.h>

namespace Util
{

   const std::string& LabelBuffer::load(std::istream& in)
   {

      if (isLoaded_) {
         UTIL_THROW("Attempt to load full LabelBuffer");
      }

      // Extract string, which may include trailing bracket
      buffer_.clear();
      in >> buffer_;
      isLoaded_ = true;

      // Copy trimmed version of buffer_ to label_
      char c = buffer_[0];
      int  n = buffer_.length(); 
      int  i = 1;
      while (i < n && c != '{') {
         label_.push_back(c);
         c = buffer_[i];
         ++i;
      }

      return label_;
      
   }
   
   const std::string& LabelBuffer::unload()
   {
      if (!isLoaded_) {
         UTIL_THROW("Attempt to unload empty LabelBuffer");
      }
      isLoaded_ = false;
      return buffer_;
   }

   bool LabelBuffer::isLoaded()
   {  return isLoaded_; }
}
