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
#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>


#include "lcio.h"
#include "IO/LCWriter.h"

#include "UTIL/LCTOOLS.h"
#include "UTIL/Operators.h"
#include "UTIL/LCObjectHandle.h"
#include "UTIL/LCTime.h"
#include "UTIL/CellIDDecoder.h"
#include "UTIL/PIDHandler.h"

#include "IMPL/LCEventImpl.h" 
#include "IMPL/LCRunHeaderImpl.h" 
#include "IMPL/LCCollectionVec.h"
#include "IMPL/SimCalorimeterHitImpl.h"
#include "IMPL/SimTrackerHitImpl.h"
#include "IMPL/MCParticleImpl.h" 
//#include "IMPL/LCFlagImpl.h" 
#include "IMPL/LCTOOLS.h"

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

void processToLegacy(podio::EventStore& store,
                     LCEventImpl* evt,
                     bool verbose,
                     int rn,
                     int eventNumber)
{

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

  auto& mcps = store.get<lcio2::MCParticleCollection>("MCParticle");
  auto& SiTrackerBarrelHits = store.get<lcio2::SimTrackerHitCollection>("SiTrackerBarrelHits");
  auto& SiTrackerEndcapHits = store.get<lcio2::SimTrackerHitCollection>("SiTrackerEndcapHits");
  auto& SiTrackerForwardHits = store.get<lcio2::SimTrackerHitCollection>("SiTrackerForwardHits");
  auto& SiVertexBarrelHits  = store.get<lcio2::SimTrackerHitCollection>("SiVertexBarrelHits");
  auto& SiVertexEndcapHits  = store.get<lcio2::SimTrackerHitCollection>("SiVertexEndcapHits");
  auto& Tracks = store.get<lcio2::TrackCollection>("Tracks");
  auto& ReconClusters = store.get<lcio2::ClusterCollection>("ReconClusters");
  auto& EM_BARREL = store.get<lcio2::CalorimeterHitCollection>("EM_BARREL");
  auto& EM_ENDCAP = store.get<lcio2::CalorimeterHitCollection>("EM_ENDCAP");
  auto& HAD_BARREL = store.get<lcio2::CalorimeterHitCollection>("HAD_BARREL");
  auto& HAD_ENDCAP = store.get<lcio2::CalorimeterHitCollection>("HAD_ENDCAP");
  auto& PandoraPFOCollection = store.get<lcio2::ReconstructedParticleCollection>("PandoraPFOCollection");
  auto& BeamCalHits = store.get<lcio2::SimCalorimeterHitCollection>("BeamCalHits");
  //auto& VXD_RawTrackerHits = store.get<lcio2::TrackerRawDataCollection>("VXD_RawTrackerHits");
  auto& HelicalTrackHits = store.get<lcio2::TrackerHitCollection>("HelicalTrackHits");


  // create and add some mc particles 
  LCCollectionVec* mcVec = new LCCollectionVec( LCIO::MCPARTICLE )  ;

  // debug only - don't write MCParticles to output file:
  // mcVec->setTransient() ;

  // debug only - add the same particle to more than one collection
  //LCCollectionVec* mcVec2 = new LCCollectionVec( LCIO::MCPARTICLE )  ;

  if( mcps.isValid() ){
    //-------- print relations for debugging:
    for(const auto& p : mcps ){

      MCParticleImpl* p_legacy = lcio2::from_lcio2(p);
      mcVec->push_back( p_legacy ) ;

      //std::cout << " particle " << p.getObjectID().index << " has daughters: " ;
      //for(auto it = p.daughters_begin(), end = p.daughters_end() ; it!=end ; ++it ){
      //  std::cout << " " << it->getObjectID().index ;
      //}
      //std::cout << "  and parents: " ;
      //for(auto it = p.parents_begin(), end = p.parents_end() ; it!=end ; ++it ){
      //  std::cout << " " << it->getObjectID().index ;
      //}
      //std::cout << std::endl ;
    }
  }

  if(verbose) {
    std::cout << "mcps: " << mcps.size() << std::endl;
  }

  //MCParticleImpl* mom = new MCParticleImpl ;
  //mom->setPDG( 1  ) ;
  //float p0[3] = { 0. , 0. , 1000. } ;
  //mom->setMomentum( p0 ) ;
  //mom->setMass( 3.01 ) ;

  //for(int j=0;j<10;j++){

    //MCParticleImpl* mcp = new MCParticleImpl ;

    //mcp->setPDG( 1000 * (j+1)  ) ;
    //float p[3] = { float(j*1.) , float(4./1024.) , float(8./1024.) } ;
    //mcp->setMomentum( p ) ;
    //mcp->setMass( .135 ) ;

    //// create and add some daughters
    //for(int k=0;k<3;k++){
    //  MCParticleImpl* d1 = new MCParticleImpl ;

    //  d1->setPDG( 1000 * (j+1) + 100 * (k+1)  ) ;
    //  float pd1[3] = {  float(k*1.) ,  float(4.1) ,  float(8.1) } ;
    //  d1->setMomentum( pd1 ) ;
    //  d1->setMass( .135 ) ;

    //  for(int l=0;l<2;l++){
    //    MCParticleImpl* d2 = new MCParticleImpl ;

    //    d2->setPDG( 1000 * (j+1) + 100 * (k+1) + 10 *  (l+1)  ) ;
    //    float pd2[3] = {  float(l*1.) ,  float(0.41) ,  float(4.1) } ;
    //    d2->setMomentum( pd2 ) ;
    //    d2->setMass( .135 ) ;

    //    double ep[3] = { 1.111111 , 2.2222222, 3.3333333 } ;
    //    d2->setEndpoint( ep ) ;
    //    //	      d2->setSimulatorStatus( 1234 ) ;
    //    d2->setCreatedInSimulation(true) ;
    //    d2->setBackscatter(true)         ;
    //    d2->setDecayedInTracker(true)    ;
    //    d2->setDecayedInCalorimeter(false);
    //    d2->setHasLeftDetector(false)     ;
    //    d2->setStopped(true)             ;

    //    d2->addParent( d1 );
    //    mcVec->push_back( d2 ) ;

    //    // debug only - add the same particle to more than one collection
    //    //mcVec2->push_back( d2 ) ;
    //  }
    //  d1->addParent( mcp );
    //  mcVec->push_back( d1 ) ;
    //}

    //mcp->addParent( mom );
    //mcVec->push_back( mcp ) ;
  //}
  //mcVec->push_back( mom ) ;

  // add all collections to the event
  evt->addCollection( mcVec , "MCParticle" ) ;

}


int main(){

  auto reader = podio::ROOTReader();
  auto store  = podio::EventStore();
  reader.openFile("example.root");
  //try {
  //  reader.openFile("example.root");
  //}
  //catch(std::runtime_error& err) {
  //  std::cerr<<err.what()<<". Quitting."<<std::endl;
  //  exit(1);
  //}
  store.setReader(&reader);

  int rn = 22;

  // ---------------------------
  // legacy-lcio output 
  LCWriter* lcWrt = LCFactory::getInstance()->createLCWriter()  ;


  // turn off compression for first run ...
  lcWrt->setCompressionLevel( 0 ) ;             
  lcWrt->open( "test.slcio" , LCIO::WRITE_NEW )  ;

  LCRunHeaderImpl* runHdr = new LCRunHeaderImpl ; 
  runHdr->setRunNumber( rn ) ;

  std::string detName("D09TileHcal")  ;
  runHdr->setDetectorName( detName ) ;

  std::stringstream description ; 
  description << " run: " << rn <<" just for testing lcio  - no physics !" ;
  runHdr->setDescription( description.str()  ) ;

  std::string ecalName("ECAL007") ;
  runHdr->addActiveSubdetector( ecalName ) ;

  std::string tpcName("TPC4711") ;
  runHdr->addActiveSubdetector( tpcName ) ;


  // add some parameters to the run header 
  //       StringVec sv1 ;
  //       sv1.push_back("simjob.cc") ;
  //       runHdr->parameters().setValues( "SimulationProgram" , sv1 ) ; 
  runHdr->parameters().setValue( "Program" , "lcio2legacy.cc" ) ; 
  IntVec iv(3) ;
  iv[0] = 1 ;
  iv[1] = 2 ;
  iv[2] = 3 ;
  runHdr->parameters().setValues( "SomeIndices" , iv ) ; 

  lcWrt->writeRunHeader( runHdr ) ;

  // event loop
  unsigned nEvents = reader.getEntries();
  for(unsigned i=0; i<nEvents; ++i) {

    bool vrb = false;
    if(i%100==0) {
      std::cout<<"reading event "<<i<<std::endl;
      vrb = true;
    }

    // we need to use the implementation classes here 
    LCEventImpl*  evt = new LCEventImpl() ;

    evt->setRunNumber(  rn   ) ;
    evt->setEventNumber( i ) ;
    LCTime now ;
    evt->setTimeStamp( now.timeStamp()  ) ;
    evt->setDetectorName( detName ) ;

    evt->setWeight( 0.5) ;

    evt->parameters().setValue("Description"," event can have it's own set of parameters" ) ;
    evt->parameters().setValue("Thrust", (float) 0.671 ) ;

    FloatVec fv ;
    fv.push_back( 1.1 ) ;
    fv.push_back( 2.2 ) ;
    fv.push_back( 3.3 ) ;
    evt->parameters().setValues( "SomeNumbers" , fv ) ; 


    processToLegacy(store, evt, vrb, rn, i);



    //  dont use this (compatibility with Fortran simjob.F)
    //  if( ! (i%100) ) cout << ". " << flush  ;

    // write the event to the file
    lcWrt->writeEvent( evt ) ;

    // dump the event to the screen 
    if(i%100==0)
      LCTOOLS::dumpEvent( evt ) ;

    // ------------ IMPORTANT ------------- !
    // we created the event so we need to delete it ...
    delete evt ;
    // -------------------------------------

    store.clear();
    reader.endOfEvent();

  } // evt loop

  delete runHdr ;

  lcWrt->close() ;
  delete lcWrt ;


  return 0;
}
