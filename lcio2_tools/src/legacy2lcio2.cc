#include "lcio2/EventInfoCollection.h"
#include "lcio2/MCParticleCollection.h"
#include "lcio2/MCParticle.h"
#include "lcio2utils.h"
#include "lcio2/TrackerHitCollection.h"
#include "lcio2/SimTrackerHitCollection.h"
#include "lcio2/CalorimeterHitCollection.h"
#include "lcio2/TrackCollection.h"
#include "lcio2/ClusterCollection.h"

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

//BeamCalHits                   SimCalorimeterHit              695 
//CalorimeterHitRelations       LCRelation                     589
//EM_BARREL                     CalorimeterHit                   4            
//EM_ENDCAP                     CalorimeterHit                 558
//HAD_BARREL                    CalorimeterHit                   0
//HAD_ENDCAP                    CalorimeterHit                  27
//HelicalTrackHitRelations      LCRelation                      47
//HelicalTrackHits              TrackerHit                      34
//HelicalTrackMCRelations       LCRelation                      34
//LumiCalHits                   SimCalorimeterHit               34
//MCInfo                        LCGenericObject                  1
//MCParameters                  LCGenericObject                  1
//MCParticle                    MCParticle                      43
//MUON_BARREL                   CalorimeterHit                   0
//MUON_ENDCAP                   CalorimeterHit                   0
//MuonBarrelHits                SimCalorimeterHit                0
//MuonEndcapHits                SimCalorimeterHit                0
//PandoraPFOCollection          ReconstructedParticle            6
//ReconClusters                 Cluster                          6
//SiTrackerBarrelHits           SimTrackerHit                   24
//SiTrackerEndcapHits           SimTrackerHit                   45
//SiTrackerForwardHits          SimTrackerHit                    5
//SiVertexBarrelHits            SimTrackerHit                    2
//SiVertexEndcapHits            SimTrackerHit                    7
//StateAtECal                   LCGenericObject                  1
//StateAtEnd                    LCGenericObject                  1
//StateAtStart                  LCGenericObject                  1
//TKR_RawTrackerHits            TrackerRawData                  89
//TKR_TrackerHits               TrackerHit                      40
//Tracks                        Track                            1
//VXD_RawTrackerHits            TrackerRawData                  26
//VXD_TrackerHits               TrackerHit                      14
//---------------------------------------------------------------------------

  auto& mcps = store.create<lcio2::MCParticleCollection>("MCParticle");
  auto& SiTrackerBarrelHits = store.create<lcio2::SimTrackerHitCollection>("SiTrackerBarrelHits");
  auto& SiTrackerEndcapHits = store.create<lcio2::SimTrackerHitCollection>("SiTrackerEndcapHits");
  auto& SiTrackerForwardHits = store.create<lcio2::SimTrackerHitCollection>("SiTrackerForwardHits");
  auto& SiVertexBarrelHits  = store.create<lcio2::SimTrackerHitCollection>("SiVertexBarrelHits");
  auto& SiVertexEndcapHits  = store.create<lcio2::SimTrackerHitCollection>("SiVertexEndcapHits");
  auto& Tracks = store.create<lcio2::TrackCollection>("Tracks");
  auto& ReconClusters = store.create<lcio2::ClusterCollection>("ReconClusters");
  auto& EM_BARREL = store.create<lcio2::CalorimeterHitCollection>("EM_BARREL");
  auto& EM_ENDCAP = store.create<lcio2::CalorimeterHitCollection>("EM_ENDCAP");
  auto& HAD_BARREL = store.create<lcio2::CalorimeterHitCollection>("HAD_BARREL");
  auto& HAD_ENDCAP = store.create<lcio2::CalorimeterHitCollection>("HAD_ENDCAP");

  writer.registerForWrite("MCParticle");
  writer.registerForWrite("SiTrackerBarrelHits");
  writer.registerForWrite("SiTrackerEndcapHits");
  writer.registerForWrite("SiTrackerForwardHits");
  writer.registerForWrite("SiVertexBarrelHits");
  writer.registerForWrite("SiVertexEndcapHits");
  writer.registerForWrite("Tracks");
  writer.registerForWrite("ReconClusters");
  writer.registerForWrite("EM_BARREL");
  writer.registerForWrite("EM_ENDCAP");
  writer.registerForWrite("HAD_BARREL");
  writer.registerForWrite("HAD_ENDCAP");

  EVENT::LCEvent* evt = 0 ;
  //evt = lcReader->readEvent(0);
  //auto collections  = evt->getCollectionNames();

  //for(const auto& v : (*collections) ) {
  //  auto a_collection = evt->getCollection(v);
  //  if( a_collection->getTypeName() == std::string("MCParticle") ) {
  //    auto& mcps = store.create<lcio2::MCParticleCollection>(v);
  //    writer.registerForWrite(v);
  //  }
  //}

  //----------- the event loop -----------
  while( (evt = lcReader->readNextEvent()) != 0  && nEvents < maxEvt ) {

    if( nEvents < 1 ){  // only dump first event
      UTIL::LCTOOLS::dumpEvent( evt ) ;
    }

    int  run_number   = evt->getRunNumber();
    int  event_number = evt->getEventNumber();
    auto collections  = evt->getCollectionNames();

    for(const auto& v : (*collections) ) {

      auto a_collection = evt->getCollection(v);

      // -----------------------------------------------------------
      // MCParticle
      if( a_collection->getTypeName() == std::string("MCParticle") ) {
        int n_elem = a_collection->getNumberOfElements();
        for(int i_elem = 0; i_elem<n_elem; i_elem++) {
          EVENT::MCParticle* elem  = dynamic_cast<EVENT::MCParticle*>( a_collection->getElementAt( i_elem ) ) ;
          auto a_mcp = lcio2::to_lcio2(elem);
          mcps.push_back(a_mcp);
        }
      }

      // SiTrackerBarrelHits
      if( v == std::string("SiTrackerBarrelHits") ) {
        int n_elem = a_collection->getNumberOfElements();
        for(int i_elem = 0; i_elem<n_elem; i_elem++) {
          EVENT::SimTrackerHit* elem  = dynamic_cast<EVENT::SimTrackerHit*>( a_collection->getElementAt( i_elem ) ) ;
          auto a_mcp = lcio2::to_lcio2(elem);
          SiTrackerBarrelHits.push_back(a_mcp);
        }
      }
      // SiTrackerEndcapHits
      if( v == std::string("SiTrackerEndcapHits") ) {
        int n_elem = a_collection->getNumberOfElements();
        for(int i_elem = 0; i_elem<n_elem; i_elem++) {
          EVENT::SimTrackerHit* elem  = dynamic_cast<EVENT::SimTrackerHit*>( a_collection->getElementAt( i_elem ) ) ;
          auto a_mcp = lcio2::to_lcio2(elem);
          SiTrackerEndcapHits.push_back(a_mcp);
        }
      }
      // SiTrackerForwardHits
      if( v == std::string("SiTrackerForwardHits") ) {
        int n_elem = a_collection->getNumberOfElements();
        for(int i_elem = 0; i_elem<n_elem; i_elem++) {
          EVENT::SimTrackerHit* elem  = dynamic_cast<EVENT::SimTrackerHit*>( a_collection->getElementAt( i_elem ) ) ;
          auto a_mcp = lcio2::to_lcio2(elem);
          SiTrackerForwardHits.push_back(a_mcp);
        }
      }
      // SiVertexBarrelHits
      if( v == std::string("SiVertexBarrelHits") ) {
        int n_elem = a_collection->getNumberOfElements();
        for(int i_elem = 0; i_elem<n_elem; i_elem++) {
          EVENT::SimTrackerHit* elem  = dynamic_cast<EVENT::SimTrackerHit*>( a_collection->getElementAt( i_elem ) ) ;
          auto a_mcp = lcio2::to_lcio2(elem);
          SiVertexBarrelHits.push_back(a_mcp);
        }
      }
      // SiVertexEndcapHits
      if( v == std::string("SiVertexEndcapHits") ) {
        int n_elem = a_collection->getNumberOfElements();
        for(int i_elem = 0; i_elem<n_elem; i_elem++) {
          EVENT::SimTrackerHit* elem  = dynamic_cast<EVENT::SimTrackerHit*>( a_collection->getElementAt( i_elem ) ) ;
          auto a_mcp = lcio2::to_lcio2(elem);
          SiVertexEndcapHits.push_back(a_mcp);
        }
      }

      //"Tracks"
      if( v == std::string("Tracks") ) {
        int n_elem = a_collection->getNumberOfElements();
        for(int i_elem = 0; i_elem<n_elem; i_elem++) {
          EVENT::Track* elem  = dynamic_cast<EVENT::Track*>( a_collection->getElementAt( i_elem ) ) ;
          auto a_mcp = lcio2::to_lcio2(elem);
          Tracks.push_back(a_mcp);
        }
      }
      //"ReconClusters");
      if( v == std::string("ReconClusters") ) {
        int n_elem = a_collection->getNumberOfElements();
        for(int i_elem = 0; i_elem<n_elem; i_elem++) {
          EVENT::Cluster* elem  = dynamic_cast<EVENT::Cluster*>( a_collection->getElementAt( i_elem ) ) ;
          auto a_mcp = lcio2::to_lcio2(elem);
          ReconClusters.push_back(a_mcp);
        }
      }
      //"EM_BARREL");
      if( v == std::string("EM_BARREL") ) {
        int n_elem = a_collection->getNumberOfElements();
        for(int i_elem = 0; i_elem<n_elem; i_elem++) {
          EVENT::CalorimeterHit* elem  = dynamic_cast<EVENT::CalorimeterHit*>( a_collection->getElementAt( i_elem ) ) ;
          auto a_mcp = lcio2::to_lcio2(elem);
          EM_BARREL.push_back(a_mcp);
        }
      }
      //"EM_ENDCAP");
      if( v == std::string("EM_ENDCAP") ) {
        int n_elem = a_collection->getNumberOfElements();
        for(int i_elem = 0; i_elem<n_elem; i_elem++) {
          EVENT::CalorimeterHit* elem  = dynamic_cast<EVENT::CalorimeterHit*>( a_collection->getElementAt( i_elem ) ) ;
          auto a_mcp = lcio2::to_lcio2(elem);
          EM_ENDCAP.push_back(a_mcp);
        }
      }
      //"HAD_BARREL");
      if( v == std::string("HAD_BARREL") ) {
        int n_elem = a_collection->getNumberOfElements();
        for(int i_elem = 0; i_elem<n_elem; i_elem++) {
          EVENT::CalorimeterHit* elem  = dynamic_cast<EVENT::CalorimeterHit*>( a_collection->getElementAt( i_elem ) ) ;
          auto a_mcp = lcio2::to_lcio2(elem);
          HAD_BARREL.push_back(a_mcp);
        }
      }
      //"HAD_ENDCAP");
      if( v == std::string("HAD_ENDCAP") ) {
        int n_elem = a_collection->getNumberOfElements();
        for(int i_elem = 0; i_elem<n_elem; i_elem++) {
          EVENT::CalorimeterHit* elem  = dynamic_cast<EVENT::CalorimeterHit*>( a_collection->getElementAt( i_elem ) ) ;
          auto a_mcp = lcio2::to_lcio2(elem);
          HAD_ENDCAP.push_back(a_mcp);
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
