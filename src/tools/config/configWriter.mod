namespace Tools
{

/**
* \defgroup Tools_ConfigWriter_Module Config File Writers
* \brief Configuration file writers
*
* The choice of a configuration file format for writing can
* be set in the command file for the MdPp postprocessing tool
* by a SET_CONFIG_WRITER command, with the syntax
* \code
*   SET_CONFIG_WRITER className
* \endcode
* where className represents the name of a subclass of the
* base class Tools::ConfigWriter, such as DdMdWriter or 
* HoomdConfigWriter.
* 
* The WRITE_CONFIG command writes a configuration file using
* the current choice of file format. The syntax for most
* (but not all) configu file writers is 
* \code
*   WRITE_CONFIG configFile
* \endcode
* where configFile is the path to the output configuration
* file, defined relative to the directory from which mdPp
* is invoked. 
*
* Some configuration file writers, such as HoomdConfigWriter,
* may require an additional auxiliary data file. The
* syntax in this case is
* is
* \code
*   WRITE_CONFIG auxiliaryFile configFile
* \endcode
* where auxiliaryFile is the path to the required auxiliary
* data file. 
*
* \ingroup Tools_Config_Module
*/

}
