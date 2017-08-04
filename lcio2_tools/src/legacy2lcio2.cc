#include "lcio2/MCParticle.h"
#include "lcio2/EventInfoCollection.h"
//#include "lcio2/JetCollection.h"
//

#include "lcio.h"
#include <stdio.h>

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

void writeEventTree(const char* FILEN = "legacy_reconstruction_datafile.slcio") {
    
  //--- create a ROOT file, a tree and a branch ...
  
  TFile* file = new TFile( "lcioEventTree.root" , "RECREATE");    
  
  TTree* tree = new TTree( "LCIO" , "lcio event data tree");
  
  std::string treeName("LCEvent") ;
  
  IMPL::LCEventImpl* treeEvt=0 ;

  std::string type("IMPL::LCEventImpl") ;
  
  TBranch* mcpBranch = tree->Branch( treeName.c_str(), 
				     type.c_str(), 
				     (void*) &treeEvt, 
				     1024, // record size 
				     199    // split level 
				     );
  

  std::cout << " loaded LCIO library and dictionary ... " << std::endl ;
  
  int nEvents = 0  ;

  int maxEvt = 10000 ; 

  IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader() ;

  lcReader->open( FILEN ) ;
  

  //----------- the event loop -----------
  EVENT::LCEvent* evt = 0 ;
  while( (evt = lcReader->readNextEvent()) != 0  && nEvents < maxEvt ) {
    
    if( nEvents < 3 )  // only dump first 3 events
      UTIL::LCTOOLS::dumpEvent( evt ) ;

    nEvents ++ ;

    treeEvt = (IMPL::LCEventImpl*) evt ;
    

    tree->Fill() ;

  }
  // -------- end of event loop -----------

  file->Write() ;
  file->Close() ;
  
  delete file ;


  std::cout << std::endl 
	    <<  "  " <<  nEvents 
	    << " events read from file: " 
	    << FILEN << std::endl  ;
  
  
  lcReader->close() ;
  
  delete lcReader ;
}

void processEvent(podio::EventStore& store, bool verbose,
                  podio::ROOTReader& reader) {

  // read event information
  //const lcio2::EventInfoCollection* evinfocoll(nullptr);
  //bool evinfo_available = store.get("EventInfo", evinfocoll);
  //if(evinfo_available) {
  //  auto evinfo = evinfocoll->at(0);

  //  if(verbose)
  //    std::cout << "event number " << evinfo.number() << std::endl;
  //}

  //// the following is commented out to test on-demand reading through Jet-Particle association,
  //// see below
  //// // read particles
  //// ParticleCollection* ptcs(nullptr);
  //// bool particles_available = store.get("GenParticle",ptcs);
  //// if (particles_available){
  ////   for(const auto& part : *ptcs) {
  ////     std::cout<<part.containerID()<<" "<<part.index()<<std::endl;
  ////   }
  //// }

  //// read jets
  //const lcio2::JetCollection* jrefs(nullptr);
  //bool jets_available = store.get("GenJet",jrefs);
  //std::vector<lcio2::ConstParticle> injets;

  //if (jets_available){
  //  if(verbose) {
  //    reader.getCollectionIDTable()->print();
  //    std::cout << "jet collection:" << std::endl;
  //  }
  //  for(const auto& jet : *jrefs){
  //    TLorentzVector lv = utils::lvFromPOD(jet.core().p4);
  //    if(verbose)
  //      std::cout << "\tjet: E=" << lv.E() << " "<<lv.Eta()<<" "<<lv.Phi()
  //                <<" npart="<<jet.particles_size()<<std::endl;
  //    for(auto part = jet.particles_begin(); part != jet.particles_end(); ++part) {
  //      if(part->isAvailable()) {
  //        if(verbose)
  //          std::cout<<"\t\tassociated "<<part->core()<<std::endl;
  //        injets.push_back(*part);
  //      }
  //    }
  //  }
  //}

  //// read particles
  //const lcio2::ParticleCollection* ptcs(nullptr);
  //bool particles_available = store.get("GenParticle", ptcs);
  //if (particles_available){
  //  std::vector<lcio2::Particle> muons;
  //  // there is probably a smarter way to get a vector from collection?

  //  if(verbose)
  //    std::cout << "particle collection:" << std::endl;
  //  for(const auto& ptc : *ptcs){
  //    if(verbose)
  //      std::cout<<"\t"<<ptc<<std::endl;
  //    if( ptc.core().pdgId == 4 ) {
  //      muons.push_back(ptc);
  //    }
  //  }
  //  // listing particles that are not used in a jet
  //  std::vector<lcio2::Particle> unused = utils::unused(*ptcs, injets);
  //  if(verbose)
  //    std::cout<<"unused particles: "<<unused.size()<<"/"<<ptcs->size()<<" "<<injets.size()<<std::endl;

  //  // computing isolation for first muon
  //  if(not muons.empty()) {
  //    const lcio2::Particle& muon = muons[0];
  //    float dRMax = 0.5;
  //    const std::vector<lcio2::Particle> incone = utils::inCone( muon.core().p4,
  //                                                        *ptcs,
  //                                                        dRMax);
  //    float sumpt = utils::sumPt(incone);
  //    if( verbose ) {
  //      std::cout<<"muon: "<<muon<<" sumpt "<<sumpt<<std::endl;
  //      std::cout<<"\tparticles in cone:"<<std::endl;
  //    }
  //    for(const auto& ptc : incone) {
  //      if( verbose )
  //        std::cout<<"\t"<<ptc<<std::endl;
  //    }
  //  }
  //}
}


int main(){
  auto reader = podio::ROOTReader();
  auto store = podio::EventStore();
  lcio2::MCParticle part;
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
  writeEventTree();
  return 0;
}
