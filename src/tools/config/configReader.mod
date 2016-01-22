namespace Tools
{

/**
* \defgroup Tools_ConfigReader_Module Config File Readers
* \brief Configuration file readers
*
* The choice of a configuration file format for reading can
* be set in the command file for the MdPp postprocessing tool
* by a SET_CONFIG_READER command, with the syntax
* \code
*   SET_CONFIG_READER className
* \endcode
* where className represents the name of a subclass of the
* base class Tools::ConfigReader, such as DdMdReader or 
* HoomdConfigReader.
* 
* The READ_CONFIG command reads a configuration file using
* the current choice of file format. The syntax for most
* (but not all) configu file readers is 
* \code
*   READ_CONFIG configFile
* \endcode
* where configFile is the path to the input configuration
* file, defined relative to the directory from which mdPp
* is invoked. 
*
* Some configuration file readers, such as HoomdConfigReader,
* may require an additional auxiliary data file. In this 
* case, they synatx of the READ_CONFIG command is
* \code
*   READ_CONFIG auxiliaryFile configFile
* \endcode
* where auxiliaryFile is the path to the required auxiliary
* data file. 
*
* \ingroup Tools_Config_Module
*/

}
