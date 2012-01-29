#ifndef LABEL_BUFFER_H
#define LABEL_BUFFER_H

#include <fstream>
#include <string>

namespace Util
{

   /**
   * A LabelBuffer stores a label to allow a one token look ahead.
   *  
   * The LabelBuffer class is used by Label and ParamComposite to
   * allow decisions to be made based on the contents of a label,
   * while leaving the label available to be read within the 
   * readParam method of the next ParamComponent in the file. 
   * The static load() method reads a string from file and stores
   * it in a std::string buffer. The unload() method returns the
   * value of the buffer, if any, or throws an Exception.
   */
   class LabelBuffer 
   {

   public:
   
      /**
      * Read a label string and store it in a buffer.
      * 
      * Throws an Exception if the buffer is already loaded.
      *
      * \return Trimmed copy of label string.
      */
      static const std::string& load(std::istream& in);

      /**
      * Return contents of buffer, and clear buffer.
      * 
      * Throws and Exception if the buffer is not loaded.
      *
      * \return untrimmed copy of buffer contents.
      */
      static const std::string& unload();

      /**
      * Does the buffer contain a string?
      */
      static bool isLoaded();

   private:

      /**
      * Buffer hold contents of string read by peek.
      */ 
      static std::string buffer_;
   
      /**
      * Trimmed label, without any trailing bracket(s).
      */ 
      static std::string label_;

      /**
      * True if the buffer contains a string.
      */
      static bool        isLoaded_;

   };
   
}
#endif
