#include "lcio2/EventInfoCollection.h"
#include "lcio2/MCParticleCollection.h"
#include "lcio2/MCParticle.h"
#include "lcio2utils.h"
#include "lcio2/TrackerHitCollection.h"
#include "lcio2/SimTrackerHitCollection.h"
#include "lcio2/CalorimeterHitCollection.h"
#include "lcio2/SimCalorimeterHitCollection.h"
#include "lcio2/TrackCollection.h"
#include "lcio2/ClusterCollection.h"
#include "lcio2/ReconstructedParticleCollection.h"
#include "lcio2/TrackerRawDataCollection.h"

#include <stdio.h>
#include <typeinfo>
#include <fstream>
#include "getopt.h"

#include "lcio.h"
#include "UTIL/LCTOOLS.h"
#include "UTIL/Operators.h"
#include "UTIL/LCObjectHandle.h"
#include "UTIL/LCTime.h"
#include "UTIL/CellIDDecoder.h"
#include "UTIL/PIDHandler.h"

#include "EVENT/LCCollection.h"
#include "EVENT/SimCalorimeterHit.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/RawCalorimeterHit.h"
#include "EVENT/SimTrackerHit.h"
#include "EVENT/TrackerHit.h"
#include "EVENT/TrackerHitPlane.h"
#include "EVENT/TrackerHitZCylinder.h"
#include "EVENT/TPCHit.h"
#include "EVENT/TrackerRawData.h"
#include "EVENT/TrackerData.h"
#include "EVENT/TrackerPulse.h"
#include "EVENT/Track.h"
#include "EVENT/LCIO.h"
#include "EVENT/MCParticle.h"
#include "EVENT/LCFloatVec.h"
#include "EVENT/LCIntVec.h"
#include "EVENT/Track.h"
#include "EVENT/Cluster.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/Vertex.h"
#include "EVENT/LCGenericObject.h"
#include "EVENT/LCRelation.h"

//#include "IMPL/LCFlagImpl.h"
#include "UTIL/BitSet32.h"
#include "UTIL/LCIOTypeInfo.h"
#include "LCIOSTLTypes.h"

#include "IO/LCReader.h"
#include "IMPL/LCTOOLS.h"
#include "EVENT/LCRunHeader.h" 

//#include "EVENT/SimCalorimeterHit.h" 
#include "EVENT/CalorimeterHit.h" 
//#include "EVENT/RawCalorimeterHit.h" 
#include "EVENT/ReconstructedParticle.h"
//#include "UTIL/CellIDDecoder.h"
#include "UTIL/Operators.h"
#include "UTIL/LCIterator.h"
// Utility functions
//#include "utilities/VectorUtils.h"
//#include "utilities/ParticleUtils.h"

// ROOT
#include "TBranch.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"

// STL
#include <vector>
#include <iostream>
#include <map>

// podio specific includes
#include "podio/EventStore.h"
#include "podio/ROOTReader.h"
#include "podio/ROOTWriter.h"

template<class T1, class T2>
std::map<std::string, T2*>  build_collection_map(LCEvent* evt, podio::EventStore& store, podio::ROOTWriter& writer)
{
  auto collections  = evt->getCollectionNames();
  std::map<std::string, T2*> res;
  for(const auto& v : (*collections) ) {
    auto a_collection = evt->getCollection(v);
    if( a_collection->getTypeName() == std::string(UTIL::lctypename<T1>()) ) {
      T2* c2 = new T2();//store.create<T2>(v);
      res[v]  = c2;
      store.registerCollection(v,c2);
      writer.registerForWrite(v);
    }
  }
  return res;
}

template<class T>
void fill_collections(LCEvent* evt, const std::string& v, T& cols_map, bool vrb = false  )
{
  if( cols_map.count(v) > 0 ){
    auto a_col = cols_map.at(v);
    auto a_collection = evt->getCollection(v);

    using col_type = decltype(a_col->create());

    int n_elem = a_collection->getNumberOfElements();

    using leg_type = decltype(lcio2::from_lcio2(std::declval<col_type>()));

    for(int i_elem = 0; i_elem<n_elem; i_elem++) {
      auto elem  = dynamic_cast<leg_type>( a_collection->getElementAt( i_elem ) ) ;
      auto a_mcp = lcio2::to_lcio2(elem);
      a_col->push_back(a_mcp);
    }
    if( vrb ){
      std::cout << std::setw(30) << v <<  " : " << a_col->size()  << std::endl;
    }
  }
}

std::string  exec(const char* cmd);
bool         fexists(const std::string& filename);
void         check_field_maps();
void         download_field_maps();
void         print_help();

int main(int argc,char** argv) {

  int          run_number        = 0;
  int          number_of_events  = -1;
  int          event_start_offset = 0;
  std::string  input_file_name   = "";//legacy_reconstruction_datafile.slcio";
  std::string  output_file_name  = "";
  std::string  output_tree_name  = "";
  std::string  output_dir        = "";
  std::string  theRest           = "";

  //---------------------------------------------------------------------------

  int index = 0;
  int iarg  = 0;
  opterr    = 1; //turn off getopt error message
  const struct option longopts[] =
  {
    {"run",         required_argument,  0, 'r'},
    {"input",       required_argument,  0, 'i'},
    {"output",      required_argument,  0, 'o'},
    {"dir",         required_argument,  0, 'D'},
    {"treename",    required_argument,  0, 't'},
    {"events",      required_argument,  0, 'n'},
    {"offset",      required_argument,  0, 'O'},
    {"help",        no_argument,        0, 'h'},
    {0,0,0,0}
  };
  while(iarg != -1) {
    iarg = getopt_long(argc, argv, "o:t:D:r:i:O:n:h", longopts, &index);

    switch (iarg)
    {
      case 'i':
        input_file_name = optarg;
        break;

      case 'r':
        run_number = atoi( optarg );
        break;

      case 'n':
        number_of_events = atoi( optarg );
        break;

      case 'O':
        event_start_offset = atoi( optarg );
        if(event_start_offset < 0 ) {
          std::cout << "Error: -O, --offset is negative! Aborting.\n";
          std::exit(EXIT_FAILURE);
        }
        break;

      case 't':
        output_tree_name = optarg;
        break;

      case 'D':
        output_dir = optarg;
        break;

      case 'o':
        output_file_name = optarg;
        if( fexists(output_file_name) ) {
          std::cout << "Error : " << output_file_name << " already exist"  << std::endl;
          exit(EXIT_FAILURE);
        }
        break;

      case 'h':
        print_help();
        exit(0);
        break;

      case '?':
        print_help();
        exit(EXIT_FAILURE);
        break;
      }
   }

   for (int i = optind; i < argc; i++) {
      theRest        += argv[i];
   }

   // Get the piped commands
   //std::vector<std::string> piped_commands;
   //if(!isatty(STDIN_FILENO)) {
   //   std::cout << "Reading piped commands...\n";
   //   std::string lineInput;
   //   while(std::getline(std::cin,lineInput)) {
   //      piped_commands.push_back(lineInput);
   //   }
   //}

  if(input_file_name == std::string("")) {
    if(optind == argc){
      print_help();
      std::cerr << "No input file supplied\n";
      exit(EXIT_FAILURE);
    }
    if(optind == argc-1) {
      input_file_name = argv[optind];
    } else {
      print_help();
      std::cerr << "Only one input file can be supplied:\n";
      std::cout << theRest << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  std::cout << input_file_name << std::endl;

  if(output_file_name == std::string("")) {
    output_file_name = input_file_name + std::string(".root");
  }


  //auto reader = podio::ROOTReader();
  auto store = podio::EventStore();
  auto writer = podio::ROOTWriter(output_file_name, &store);

  int nEvents = 0  ;
  int maxEvt = 10000 ; 
  IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader() ;
  lcReader->open( input_file_name ) ;


  // -----------------------------------------------
  // Get the first event to build collections
  EVENT::LCEvent* evt = lcReader->readNextEvent() ;
  run_number   = evt->getRunNumber();
  int  event_number = evt->getEventNumber();

  auto mcp_map                   = build_collection_map<EVENT::MCParticle, lcio2::MCParticleCollection>(evt, store, writer);
  auto SimTrackerHit_map         = build_collection_map<EVENT::SimTrackerHit, lcio2::SimTrackerHitCollection>(evt, store, writer);
  auto TrackerHit_map            = build_collection_map<EVENT::TrackerHit, lcio2::TrackerHitCollection>(evt, store, writer);
  auto SimCalorimeterHit_map     = build_collection_map<EVENT::SimCalorimeterHit, lcio2::SimCalorimeterHitCollection>(evt, store, writer);
  auto CalorimeterHit_map        = build_collection_map<EVENT::CalorimeterHit, lcio2::CalorimeterHitCollection>(evt, store, writer);
  auto Track_map                 = build_collection_map<EVENT::Track, lcio2::TrackCollection>(evt, store, writer);
  auto Cluster_map               = build_collection_map<EVENT::Cluster, lcio2::ClusterCollection>(evt, store, writer);
  auto ReconstructedParticle_map = build_collection_map<EVENT::ReconstructedParticle, lcio2::ReconstructedParticleCollection>(evt, store, writer);

  // close and open so the first event is read again
  lcReader->close();
  lcReader->open( input_file_name ) ;

  //----------- the event loop -----------
  while( (evt = lcReader->readNextEvent()) != 0  && nEvents < maxEvt ) {

    bool vrb = false;
    if(nEvents%100==0) {
      std::cout<<"reading event "<<nEvents<<std::endl;
      vrb = true;
    }
    if( vrb ){
      UTIL::LCTOOLS::dumpEvent( evt ) ;
    }

    run_number   = evt->getRunNumber();
    event_number = evt->getEventNumber();
    auto collections  = evt->getCollectionNames();

    for(const auto& v : (*collections) ) {

      //auto a_collection = evt->getCollection(v);

      fill_collections( evt, v, mcp_map, vrb );
      fill_collections( evt, v, SimTrackerHit_map, vrb );
      fill_collections( evt, v, TrackerHit_map           , vrb );
      fill_collections( evt, v, SimCalorimeterHit_map    , vrb );
      fill_collections( evt, v, CalorimeterHit_map       , vrb );
      fill_collections( evt, v, Track_map                , vrb );
      fill_collections( evt, v, Cluster_map              , vrb );
      fill_collections( evt, v, ReconstructedParticle_map, vrb );


    }

    nEvents ++ ;

    writer.writeEvent();
    store.clearCollections();
  }
  // -------- end of event loop -----------


  std::cout << std::endl 
    <<  "  " <<  nEvents 
    << " events read from file: " 
    << input_file_name << std::endl  ;

  lcReader->close() ;
  delete lcReader ;

  writer.finish();

  return 0;
}

//______________________________________________________________________________

std::string exec(const char* cmd) {
   std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
   if (!pipe) return "ERROR";
   char buffer[128];
   std::string result = "";
   while (!feof(pipe.get())) {
      if (fgets(buffer, 128, pipe.get()) != NULL)
         result += buffer;
   }
   return result;
}
//______________________________________________________________________________

bool fexists(const std::string& filename) {
   std::ifstream ifile(filename.c_str());
   if( ifile ) return true;
   return false;
}
//______________________________________________________________________________


void check_field_maps()
{
   //wget http://clasweb.jlab.org/12gev/field_maps/clas12SolenoidFieldMap.dat
   //wget http://clasweb.jlab.org/12gev/field_maps/clas12TorusOriginalMap.dat
   //bool failed = false;

   //if( ! fexists(C12SIM_DATA_DIR"/clas12SolenoidFieldMap.dat") ) {
   //   std::cerr << "Error: Solenoid field map missing!" << std::endl;
   //   std::cerr << "file \"" C12SIM_DATA_DIR "/clas12SolenoidFieldMap.dat\" not found. " << std::endl;
   //   failed = true;
   //}
   //if( ! fexists(C12SIM_DATA_DIR"/clas12TorusOriginalMap.dat") ) {
   //   std::cerr << "Error: Torus field map missing!" << std::endl;
   //   std::cerr << "file \"" C12SIM_DATA_DIR "/clas12TorusOriginalMap.dat\" not found. " << std::endl;
   //   failed = true;
   //}
   //if(failed) {
   //   std::cerr << "Use \"c12sim --dl-field-maps\" to download and install the field maps" << std::endl;
   //   exit(0);
   //}
}

void download_field_maps()
{
   //wget http://clasweb.jlab.org/12gev/field_maps/clas12SolenoidFieldMap.dat
   //wget http://clasweb.jlab.org/12gev/field_maps/clas12TorusOriginalMap.dat

   //std::cout << "c12sim field map directory : " C12SIM_DATA_DIR << std::endl;
   //if( ! fexists(C12SIM_DATA_DIR "/clas12SolenoidFieldMap.dat") ) {
   //   std::cout << " Downloading  http://clasweb.jlab.org/12gev/field_maps/clas12SolenoidFieldMap.dat" << std::endl;
   //   exec(" cd " C12SIM_DATA_DIR " ; wget http://clasweb.jlab.org/12gev/field_maps/clas12SolenoidFieldMap.dat && pwd && ls -lrth ");
   //}
   //if( ! fexists(C12SIM_DATA_DIR "/clas12TorusOriginalMap.dat") ) {
   //   std::cout << " Downloading  http://clasweb.jlab.org/12gev/field_maps/clas12TorusOriginalMap.dat" << std::endl;
   //   exec(" cd " C12SIM_DATA_DIR " ; wget http://clasweb.jlab.org/12gev/field_maps/clas12TorusOriginalMap.dat && pwd && ls -lrth ");
   //}
}
//______________________________________________________________________________

void print_help() {

   std::cout << "usage: lcio2legacy [options] file.root    \n";
   std::cout << "       Converts podio based lcio2 data into legacy-lcio format.\n";
   std::cout << "Options:                               \n";
   std::cout << "    -r, --run=NUMBER         Set simulation \"run\" NUMBER\n";
   std::cout << "    -i, --input=FILE         Set the input file from which events will be read\n";
   std::cout << "    -o, --output=NAME        Set the output file name which will have the run number appended\n";
   std::cout << "                             Default is \"clas12sim\". Note this is just the file basename; use -D to set the directory \n";
   std::cout << "    -n, --events=#           Causes the execution of the ui command \"/run/beamOn #\".\n";
   std::cout << "    -O, --offset=#           Starting event offset.\n";
   std::cout << "    -D, --dir=NAME           Set the output directory. The default is \"data/rootfiles/\"\n";
   std::cout << "    -t, --treename=NAME      Set the output tree name. The default is \"clasdigi_hits\"\n";
   //   {"solenoid-field", required_argument,  0, 'S'},
   //   {"toroid-field", required_argument,  0, 'T'},
   //   {"torus-field", required_argument,  0, 'T'},
}
//______________________________________________________________________________

