#include <unistd.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <memory>
#include <vector>
#include <algorithm>

#include "xmlParser.h"
#include <stdio.h>


/**
* Tool to extract topology information from MOL2 or Hoomd XML file format
* and use it to generate a Simpatico parameter file
*/

using namespace std;

class Bond {
  public:
 
  Bond() {};

  Bond(unsigned long atomAIn, unsigned long atomBIn)
  {
    atomA = atomAIn;
    atomB = atomBIn;
  }
                                                   
  unsigned long atomA, atomB;

  bool operator == (const Bond& other) const
  {
    if ((this->atomA == other.atomA) && (this->atomB == other.atomB))
      return true;
    else return false;
  }

//  Bond & operator = (const Bond& other) 
//  {
//    if (this != &other)
//    {
//      atomA = other.atomA;
//      atomB = other.atomB;
//    }
//    return *this;
//  }
};

class Molecule 
   {
   public: 

   Molecule() {
     atomTypeIds = new vector<int>;
     bonds = new vector<Bond>;
   };

   ~Molecule() {
     delete(atomTypeIds);
     delete(bonds);
   }
   vector<int> *atomTypeIds;
   vector<Bond> *bonds;

   bool operator == (const Molecule& other) const
   {
     if ((*(this->bonds) == *(other.bonds))
       && (*(this->atomTypeIds) == *(other.atomTypeIds))) return true;
     else return false;
   }

   };

class Species
   {
   public:
   Species() : nMolecule(0) {
   };

   Molecule *molecule;
   unsigned long nMolecule;
   };

class TopologyFile
   {
   public:
     TopologyFile() : uniqueAtomTypes(0), lx(0.0), ly(0.0), lz(0.0) {}

     vector<string> atomTypes;
     vector<int> atomTypeIds;
     vector<string> atomTypeNames;
     vector<Bond> bonds;
     int uniqueAtomTypes;
     double lx,ly,lz; // Box dimensions
   };


// This function reads a HOOMD XML file
// portions of the code are taken from the the HOOMDInitializer class (HOOMD source code)
int ReadHOOMDFile(TopologyFile &topologyFile,string filename)
{
   string line;

   cout << "Reading " <<filename << "..." << endl;
   XMLResults results;
   XMLNode rootNode = XMLNode::parseFile(filename.c_str(),"hoomd_xml", &results);

   if (results.error != eXMLErrorNone)
   {
     if (results.error==eXMLErrorFirstTagNotFound)
     {
        cout << "Error reading XML file! Root node of " << filename << " is not <hoomd_xml>" << endl;
     } else {
      ostringstream errorMessage;
      errorMessage << XMLNode::getError(results.error) << " in file "
      << filename << " at line " << results.nLine << " col "
      << results.nColumn;
      cout << "Error reading XML file: " << errorMessage.str() << endl;
     }
     return 1;
   }

  string xmlVersion;
  if (rootNode.isAttributeSet("version"))
     xmlVersion = rootNode.getAttribute("version");

  vector<string> valid_versions;
  valid_versions.push_back("1.3");
  bool valid = false;
  vector<string>::iterator i;
  for (i = valid_versions.begin(); i != valid_versions.end(); ++i)
  {
    if (xmlVersion == *i)
    {
      valid = true;
      break;
    }
  }
  if (!valid)
    cout << "Warning: hoomd_xml version not 1.3! Continuing anyway." << endl;
 
  int numConfigurations = rootNode.nChildNode("configuration");  
  if (numConfigurations != 1) {
    cout << "Error: Invalid number (" << numConfigurations << ") of configurations found!"
         << endl;
    return 1;
  }

  // extract configuration node
  XMLNode configurationNode = rootNode.getChildNode("configuration");

  // loop through all child nodes of the configuration
  for (int cur_node=0; cur_node < configurationNode.nChildNode(); cur_node++)
  {
    // extract the name
    XMLNode node = configurationNode.getChildNode(cur_node);
    string name = node.getName();
    transform(name.begin(), name.end(), name.begin(), ::tolower);

    if (name == "box") {
      // Box dimensions
      double lx,ly,lz;

      istringstream iss;
      if (!node.isAttributeSet("lx"))
      {
        cerr << "Error: lx not set in <box> node" << endl;
        return 1;
      }
      iss.str(node.getAttribute("lx"));
      iss >> lx;
      iss.clear();

      if (!node.isAttributeSet("ly"))
      {
        cerr << "Error: ly not set in <box> node" << endl;
        return 1;
      }
      iss.str(node.getAttribute("ly"));
      iss >> ly;
      iss.clear();

      if (!node.isAttributeSet("lz"))
      {
        cerr << "Error: lx not set in <box> node" << endl;
        return 1;
      }
      iss.str(node.getAttribute("lz"));
      iss >> lz;
      iss.clear();

      topologyFile.lx = lx;
      topologyFile.ly = ly;
      topologyFile.lz = lz;
    } else if (name == "bond") {
      // Bond list
      string allText;

      for (int i = 0; i < node.nText(); i++)
         allText += string(node.getText(i)) + string ("\n");    

      istringstream parser;
      parser.str(allText);

      while (parser.good())
      {
        string typeName;
        unsigned long atomA,atomB;

        parser >> typeName >> atomA >> atomB;
        if (parser.good())
          topologyFile.bonds.push_back(Bond(atomA,atomB));
      }
    } else if (name == "type") {
      // Type list
      
      string allText;
      for (int i = 0; i < node.nText(); i++)
         allText += string(node.getText(i)) + string ("\n");

      stringstream parser;
      parser.str(allText);

      while (parser.good())
      {
        string typeName;

        parser >> typeName;
        if (parser.good())
          topologyFile.atomTypes.push_back(typeName);
      }
    }
    // skip unneeded XML nodes
  } 

  // Currently we require bonds between particles
  if (!topologyFile.bonds.size()) {
    cout << "Error: no bond information found in XML file." << endl;
    return 1;
  }

  // We require type information
  if (!topologyFile.atomTypes.size()) {
    cout << "Error: no atom type information found in XML file. Need to enable 'type' export option in Hoomd." << endl;
    return 1;
  }

  return 0;
}

int ReadMOL2File(TopologyFile &topologyFile,string filename)
{
   string line;
   unsigned long nAtom,nBond;

   // Open topology file for reading
   ifstream file;
   file.open(filename.c_str());
   if (file.fail()) {
      std::string message = "Error opening topology file. Filename: ";
      message += filename;
      std::cout << message << std::endl;
      return 1;
   }

   // Find "Molecule" RTI (Record Type Information)
   bool foundMoleculeRTI = false;
   while (file.good() && !foundMoleculeRTI )
   {
      getline(file,line);
      if (line == "@<TRIPOS>MOLECULE") {
        getline(file, line); // mol_name
        getline(file, line);
        stringstream s(line);
        s >> nAtom;
        s >> nBond;
        foundMoleculeRTI = true;
      }
   }

   if (!foundMoleculeRTI) {
     cout << "Error: did not find MOLECULE record!" << endl;
     return 1;
   }

   // resize vectors
   topologyFile.atomTypes.resize(nAtom);
   topologyFile.bonds.resize(nBond);

   // rewind file
   file.seekg(0); 

   // Read atom types
   bool foundAtomRTI = false;
   while (file.good() && !foundAtomRTI )
   {
      getline(file,line);
      if (line == "@<TRIPOS>ATOM") {
        for (unsigned long i=0; i<nAtom; i++) {
          getline(file, line);

          unsigned long atomId;
          string atomName;
          string atomType;
          float x,y,z;

          stringstream s(line);
          s >> atomId;
          s >> atomName;
          s >> x;
          s >> y;
          s >> z;
          s >> atomType;

          topologyFile.atomTypes[i]=atomType;
        }
        foundAtomRTI = true;
      }
   }


   // Read bonds
   bool foundBondRTI = false;

   // rewind file
   file.seekg(0); 

   while (file.good() && !foundBondRTI )
   {
      getline(file,line);
      if (line == "@<TRIPOS>BOND") {
        foundBondRTI=true;
        unsigned long i =0;

        while (i<nBond) {

          getline(file, line);

          unsigned long bondId;
          unsigned long atomA;
          unsigned long atomB;
          string bondType;

          stringstream s(line);
          s >> bondId;
          s >> atomA;
          s >> atomB;
          s >> bondType;

          topologyFile.bonds[i]=Bond(atomA,atomB);
          i++;
        }
      }
   }

   // Close file
   file.close();
   return 0;
}

int main(int argc, char **argv)
{

   if (argc<3) {
      cout << "Usage:" << endl;
      cout << "extract_topology [-T <inputfileformat>] <infile> <outfile>" 
           <<endl; 
      return 1;
   }

   enum fileFormats {HOOMD,MOL2};
   fileFormats fileFormat = HOOMD;
   string fileFormatName;

   opterr = 0;
   char c;
   string::iterator it;

   while ((c = getopt(argc,argv,"T:")) != -1) {
      switch (c) {
      case 'T':
        fileFormatName = string(optarg);

        it = fileFormatName.begin();
        while (it != fileFormatName.end()) {
          *it = toupper((unsigned char)(*it));
          ++it;
         }

        if (fileFormatName == "MOL2") 
          fileFormat = MOL2;
        else if (fileFormatName == "HOOMD")
          fileFormat = HOOMD;
        else {
          cout << "Invalid file format!" << endl;
          return 1;
        }
        break;
      case '?':
        std::cout << "Unknown option -" << optopt << std::endl;
        return 1;
      }
   }

   string inputFilename;
   string outputFilename;

   inputFilename = string(argv[argc-2]);
   outputFilename = string(argv[argc-1]);

   TopologyFile topologyFile;

   int ret = 0;
   switch (fileFormat)
   {
     case MOL2:
       ret = ReadMOL2File(topologyFile,inputFilename);
       break;
     case HOOMD:
       ret = ReadHOOMDFile(topologyFile,inputFilename);
       break;
     default:
       // It shouldn't go here
       return 1;
    }
   if (ret) {
     cout << "Aborting." << endl;
     return 1;
   }

   topologyFile.atomTypeIds.resize(topologyFile.atomTypes.size());
   for (unsigned long i=0; i<topologyFile.atomTypes.size(); i++) {
     bool found=false;
     int atomTypeId =0;

     string atomType = topologyFile.atomTypes[i];

     for (unsigned long int j=0; j<i && !found; j++) {
       if (topologyFile.atomTypes[j] == atomType) {
         found=true;
         atomTypeId=topologyFile.atomTypeIds[j];
       }
     }

     if (!found) {
       atomTypeId = topologyFile.uniqueAtomTypes++;
       topologyFile.atomTypeNames.push_back(atomType);
     }

     topologyFile.atomTypeIds[i]=atomTypeId;
   }


   // Read bonds and extract species information
   vector<Species *> species;
   Molecule *thisMolecule= new Molecule;
   Species *thisSpecies;

   unsigned long startAtomId=0;
   bool readFirstBond = false;
   unsigned long firstAtomId=0;
   unsigned long lastAtomB=0;

   unsigned long iBond=0;
   while (iBond<topologyFile.bonds.size()) {

     unsigned long atomA, atomB;
     atomA=topologyFile.bonds[iBond].atomA;
     atomB=topologyFile.bonds[iBond].atomB;

     if (!readFirstBond) {
        startAtomId = atomA; 
        firstAtomId = atomA;
        readFirstBond=true;
     }
     if (atomB != atomA+1) {
       cout << "Non-consecutive bonds found! "
            << "(BondId " << iBond << " between atoms " << atomA
            << " and " << atomB << ")" << endl;
       return 1;
     }

     if (thisMolecule->bonds->size()) {
       if (atomA != lastAtomB) {
         // New molecule

         // complete old molecule
         thisMolecule->atomTypeIds->push_back(topologyFile.atomTypeIds[lastAtomB-firstAtomId]);

         // Currently, all atoms need to have at least one bond
         if (atomA != lastAtomB +1) {
           cout << "Non-consecutive molecules found! "
                << "(at BondId " << iBond << ")" << endl;
           return 1;
         }
                

         // Compare thisMolecule to current species
         if (!species.size() || !((*thisMolecule)==*(thisSpecies->molecule))) {
           // New species
           thisSpecies = new Species;
           cout << "Found new species " << 
                distance(species.begin(),species.end())+1
                << " (linear "
                << thisMolecule->atomTypeIds->size() << ")" << endl;
           thisSpecies->nMolecule++;
           thisSpecies->molecule = thisMolecule;
           species.push_back(thisSpecies);
         } else {
           // Old species
           (*species.rbegin())->nMolecule++;
         }

         startAtomId = atomA;
         thisMolecule = new Molecule;
       }
     }  

    lastAtomB = atomB;
    Bond thisBond;
    thisBond.atomA = atomA-startAtomId;
    thisBond.atomB = atomB-startAtomId;
    thisMolecule->bonds->push_back(thisBond);
    thisMolecule->atomTypeIds->push_back(topologyFile.atomTypeIds[atomA-firstAtomId]);
    iBond++;
  }

  // Complete species
  (*species.rbegin())->nMolecule++;

   cout << "Found " << species.size() << " distinct species" <<
        " and " << topologyFile.uniqueAtomTypes << " unique atom types" << endl;

   // Open output file for writing
   ofstream *outputFilePtr = new std::ofstream();
   outputFilePtr->open(outputFilename.c_str());
   if (outputFilePtr->fail()) {
      std::string message = "Error opening file for writing. Filename: ";
      message += outputFilename;
      std::cout << message << std::endl;
      return 1;
   }

   vector<Species *>::const_iterator i;

   *outputFilePtr << "MdSimulation{" << endl;

   *outputFilePtr << "  FileMaster{" << endl 
                  << "    commandFileName commands" << endl
                  << "    inputPrefix ./" << endl
                  << "    outputPrefix ./" << endl
                  << "  }" << endl;
   *outputFilePtr << "  nAtomType " << topologyFile.uniqueAtomTypes << endl
                  << "  nBondType " << 1 << endl;
   
   *outputFilePtr << "  atomTypes "; 
   for (int j=0; j<topologyFile.uniqueAtomTypes; j++)
     *outputFilePtr << topologyFile.atomTypeNames[j] << " 1.0" << endl;

   *outputFilePtr << "  maskedPairPolicy MaskBonded" << endl
                  << endl
                  << "  SpeciesManager{" << endl;

   for (i=species.begin(); i!=species.end(); ++i)
   {   
     *outputFilePtr << "    Species{" << endl;
     *outputFilePtr << "      moleculeCapacity "
                    << (*i)->nMolecule << endl;
     *outputFilePtr << "      nAtom " 
                    << (*i)->molecule->atomTypeIds->size() << endl;
     *outputFilePtr << "      nBond "
                    << (*i)->molecule->bonds->size() << endl;
     *outputFilePtr << "      atomTypeIds";
     vector<int>::const_iterator j;
     for (j=(*i)->molecule->atomTypeIds->begin(); j!=(*i)->molecule->atomTypeIds->end(); ++j)
     {
       *outputFilePtr << " " << *j;
     }
     *outputFilePtr << endl;
     *outputFilePtr << "      speciesBonds"; 

     vector<Bond>::const_iterator b;
     for (b=(*i)->molecule->bonds->begin(); b!=(*i)->molecule->bonds->end(); ++b) {
       *outputFilePtr << " " << b->atomA << " " << b->atomB << " 0" << endl;
     }
     *outputFilePtr << "    }" << endl;
   }
   *outputFilePtr << "  }" << endl
                  << "  Random{" << endl
                  << "    seed 0" << endl
                  << "  }" << endl;

   // Initialize the parameter file with some standard values
   // these have to be adjusted according to the actual parameters
   // before running certain diagnostics
   *outputFilePtr << "  MdSystem{" << endl;
   *outputFilePtr << "    EnergyEnsemble{" << endl;
   *outputFilePtr << "      type                  isothermal" << endl;
   *outputFilePtr << "      temperature           1.000000000000e+00" << endl;
   *outputFilePtr << "    }" << endl;
   *outputFilePtr << "    BoundaryEnsemble{" << endl;
   *outputFilePtr << "      type                  rigid" << endl;
   *outputFilePtr << "    }" << endl;
   *outputFilePtr << "    maxBoundary             orthorhombic "
                  << topologyFile.lx << " "
                  << topologyFile.ly << " "
                  << topologyFile.lz << endl;

   *outputFilePtr << "    LJPair{" << endl;

   // 2D Matrix of epsilon values
   *outputFilePtr << "      epsilon" <<endl;
   for (int j=0; j<topologyFile.uniqueAtomTypes; j++) {
     for (int k=0; k< topologyFile.uniqueAtomTypes; k++)
       *outputFilePtr << "1.00000 ";
     *outputFilePtr << endl;
   }

   // 2D Matrix of sigma values
   *outputFilePtr << "      sigma" << endl;
   for (int j=0; j<topologyFile.uniqueAtomTypes; j++) {
     for (int k=0; k< topologyFile.uniqueAtomTypes; k++)
       *outputFilePtr << "1.00000 ";
     *outputFilePtr << endl;
   }  

   // 2D Matrix of cutoff values
   *outputFilePtr << "      cutoff" << endl;
   for (int j=0; j<topologyFile.uniqueAtomTypes; j++) {
     for (int k=0; k< topologyFile.uniqueAtomTypes; k++)
       *outputFilePtr << "1.12246 ";
     *outputFilePtr << endl;
   }
   *outputFilePtr << "    }" << endl;

   *outputFilePtr << "    HarmonicBond{" << endl;
   *outputFilePtr << "      kappa                 4.000000000000e+02" << endl;
   *outputFilePtr << "      length                1.000000000000e+00" << endl;
   *outputFilePtr << "    }" << endl;
   *outputFilePtr << "    PairList{" << endl;
   *outputFilePtr << "      atomCapacity          ";

   *outputFilePtr << topologyFile.atomTypes.size() << endl;
   //  standard value: ten neighbors on average 
   *outputFilePtr << "      pairCapacity          " << topologyFile.atomTypes.size()*10 << endl;
   *outputFilePtr << "      skin                  3.000000000000e-01" << endl;
   *outputFilePtr << "    }" << endl;
   *outputFilePtr << "    NVTIntegrator{" <<endl;
   *outputFilePtr << "      dt                    1.200000000000e-02" << endl;
   *outputFilePtr << "      tauT                  1.000000000000" << endl;
   *outputFilePtr << "    }" << endl;
   *outputFilePtr << "  }" << endl;
   *outputFilePtr << "  DiagnosticManager{" << endl;
   *outputFilePtr << "    baseInterval                         1" << endl;
   *outputFilePtr << endl;
   *outputFilePtr << "    LogProgress{" << endl;
   *outputFilePtr << "      interval                          1000" << endl;
   *outputFilePtr << "    }" << endl;
   *outputFilePtr << "  }" << endl;
   *outputFilePtr << "}" << endl;

   // Normal completion
   return 0;

}
