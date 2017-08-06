#include "lcio2/EventInfoCollection.h"
#include "lcio2/MCParticleCollection.h"
#include "lcio2/MCParticle.h"
#include "lcio2utils.h"

#include "lcio.h"
#include <stdio.h>
#include <typeinfo>

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

// podio specific includes
#include "podio/EventStore.h"
#include "podio/ROOTReader.h"
#include "podio/ROOTWriter.h"

//template<class T>
//void setMomentum(T& p, double pt, double eta, double phi) {
//    p.Px = TMath::Cos(phi)*pt;
//    p.Py = TMath::Sin(phi)*pt;
//    p.Pz = TMath::Sinh(eta)*pt;
//}
//// First get the collection (note that this will work differently in FCCSW, see below)
//auto& myParticleColl = store.create<fcc::ParticleCollection>("recoParticles");
//auto& myMCParticleCollection = store.create<fcc::MCParticleCollection>("genParticles");
//auto particle = myParticleCollection.create();
//auto mcParticle = myMCParticleCollection.create();
//setKinematics(particle.Core, 20, 0.1, 0.5);
//setKinematics(mcParticle.Core, 20, 0.1, 0.5);

void writeEventTree(podio::EventStore& store,
                    bool verbose,
                    podio::ROOTWriter& writer,
                    char* FILEN = "legacy_reconstruction_datafile.slcio")
{

  //--- create a ROOT file, a tree and a branch ...

  //TFile* file = new TFile( "lcioEventTree.root" , "RECREATE");    
  //TTree* tree = new TTree( "LCIO" , "lcio event data tree");
  //std::string treeName("LCEvent") ;
  IMPL::LCEventImpl* treeEvt=0 ;
  //std::string type("IMPL::LCEventImpl") ;
  //TBranch* mcpBranch = tree->Branch( treeName.c_str(), 
  //      			     type.c_str(), 
  //      			     (void*) &treeEvt, 
  //      			     1024, // record size 
  //      			     199    // split level 
  //      			     );
  //
  //std::cout << " loaded LCIO library and dictionary ... " << std::endl ;

  int nEvents = 0  ;
  int maxEvt = 10000 ; 
  IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader() ;
  lcReader->open( FILEN ) ;

  auto& mcps = store.create<lcio2::MCParticleCollection>("MCParticle");
  writer.registerForWrite("MCParticle");

  //----------- the event loop -----------
  EVENT::LCEvent* evt = 0 ;
  while( (evt = lcReader->readNextEvent()) != 0  && nEvents < maxEvt ) {

    if( nEvents < 1 ){  // only dump first event
      UTIL::LCTOOLS::dumpEvent( evt ) ;
    }

    int  run_number   = evt->getRunNumber();
    int  event_number = evt->getEventNumber();
    auto collections  = evt->getCollectionNames();

    for(const auto& v : (*collections) ) {

      // MCParticle
      auto a_collection = evt->getCollection(v);
      if( a_collection->getTypeName() == std::string("MCParticle") ) {

        int n_elem = a_collection->getNumberOfElements();
        for(int i_elem = 0; i_elem<n_elem; i_elem++) {

          auto part  = dynamic_cast<MCParticle*>( a_collection->getElementAt( i_elem ) ) ;
          auto a_mcp = lcio2::to_lcio2(part);
          
          //auto a_mcp = lcio2::MCParticle();
          //a_mcp.momentum({part->getMomentum()[0] ,  part->getMomentum()[1] ,   part->getMomentum()[2] });
          //a_mcp.vertex({part->getVertex()[0] ,  part->getVertex()[1] ,   part->getVertex()[2] });
          //a_mcp.pdg( part->getPDG() );
          //a_mcp.genstatus( part->getGeneratorStatus() );
          //a_mcp.mass( part->getMass() );
          //a_mcp.charge( part->getCharge() );
          //a_mcp.time( part->getTime() );
          //a_mcp.endpoint({part->getEndpoint()[0] ,  part->getEndpoint()[1] ,   part->getEndpoint()[2] });
          //mcps.push_back(a_mcp);
          
          std::cout << a_mcp << std::endl;
        }
      }
    }

    nEvents ++ ;

    treeEvt = (IMPL::LCEventImpl*) evt ;

    //tree->Fill() ;
    writer.writeEvent();
    store.clearCollections();
  }
  // -------- end of event loop -----------

  //file->Write() ;
  //file->Close() ;
  //delete file ;

  std::cout << std::endl 
    <<  "  " <<  nEvents 
    << " events read from file: " 
    << FILEN << std::endl  ;

  lcReader->close() ;
  delete lcReader ;
}

void processEvent(podio::EventStore& store, bool verbose,
                  podio::ROOTReader& reader) {

}


int main(){

  //auto reader = podio::ROOTReader();
  auto store = podio::EventStore();
  auto writer = podio::ROOTWriter("example.root", &store);
  //lcio2::MCParticle part;
  //try {
  //  reader.openFile("example.root");
  //}
  //catch(std::runtime_error& err) {
  //  std::cerr<<err.what()<<". Quitting."<<std::endl;
  //  exit(1);
  //}
  //store.setReader(&reader);

  //bool verbose = true;

  //// unsigned nEvents = 5;
  //unsigned nEvents = reader.getEntries();
  //for(unsigned i=0; i<nEvents; ++i) {
  //  if(i%1000==0) {
  //    std::cout<<"reading event "<<i<<std::endl;
  //  }
  //  if(i>10) {
  //    verbose = false;
  //  }
  //  processEvent(store, verbose, reader);
  //  store.clear();
  //  reader.endOfEvent();
  //}
  //
  writeEventTree(store, 0, writer);


writer.finish();
  return 0;
}
