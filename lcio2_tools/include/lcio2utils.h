#ifndef LCIO2TOOLS_HH
#define LCIO2TOOLS_HH

#include "EVENT/MCParticle.h"
#include "EVENT/Cluster.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/SimCalorimeterHit.h"
#include "EVENT/RawCalorimeterHit.h"
#include "EVENT/TrackerHit.h"
#include "EVENT/SimTrackerHit.h"
//#include "EVENT/RawTrackerHit.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/Vertex.h"

#include "lcio2/MCParticleCollection.h"
#include "lcio2/MCParticle.h"
#include "lcio2/Cluster.h"
#include "lcio2/CalorimeterHit.h"
#include "lcio2/SimCalorimeterHit.h"
#include "lcio2/RawCalorimeterHit.h"
#include "lcio2/TrackerHit.h"
#include "lcio2/SimTrackerHit.h"
//#include "lcio2/RawTrackerHit.h"
#include "lcio2/ReconstructedParticle.h"
#include "lcio2/Vertex.h"


namespace lcio2 {
lcio2::MCParticle to_lcio2(EVENT::MCParticle* mcp)
  {
    return lcio2::MCParticle(
      mcp->getPDG(),
      mcp->getGeneratorStatus(),
      {mcp->getVertex()[0],mcp->getVertex()[1],mcp->getVertex()[2]},
      mcp->getCharge(),
      mcp->getMass(),
      mcp->getTime(),
      {mcp->getEndpoint()[0],mcp->getEndpoint()[1],mcp->getEndpoint()[2]},
      mcp->vertexIsNotEndpointOfParent(),
      {mcp->getMomentum()[0],mcp->getMomentum()[1],mcp->getMomentum()[2]});
  }

  Cluster to_lcio2(EVENT::Cluster* c)
  {
    //Cluster(float energy,std::array<float, 3> position,std::array<float, 6> positionError,
    //        float theta,float phi,std::array<float, 3> directionError);
    /// TODO: Fix data set to zero.
    return Cluster(
      c->getEnergy(),
      {c->getPosition()[0],c->getPosition()[1],c->getPosition()[2]},
      {0,0,0, 0,0,0},
      c->getITheta(),
      c->getIPhi(),
      {0,0,0} );
  }

  SimCalorimeterHit to_lcio2(EVENT::SimCalorimeterHit* c)
  {
    return SimCalorimeterHit(
      c->getCellID0(),
      c->getCellID1(),
      c->getEnergy(),
      c->getTimeCont(0), // FIXME:  how to treat properly time contributions?
      {c->getPosition()[0],c->getPosition()[1],c->getPosition()[2]} );
  }

  RawCalorimeterHit to_lcio2(EVENT::RawCalorimeterHit* c)
  {
    return RawCalorimeterHit(
      c->getCellID0(),
      c->getCellID1(),
      c->getAmplitude(),
      c->getTimeStamp());
  }

  CalorimeterHit to_lcio2(EVENT::CalorimeterHit* c)
  {
    return CalorimeterHit(
      c->getCellID0(),
      c->getCellID1(),
      c->getEnergy(),
      c->getTime(),
      {c->getPosition()[0],c->getPosition()[1],c->getPosition()[2]},
      c->getType());
  }

  SimTrackerHit to_lcio2(EVENT::SimTrackerHit* h)
  {
    return SimTrackerHit(
      h->getCellID0(),
      h->getCellID1(),
      {h->getPosition()[0],h->getPosition()[1],h->getPosition()[2]},
      h->getTime(),
      h->getPathLength(),
      h->getEDep(),
      {h->getMomentum()[0],h->getMomentum()[1],h->getMomentum()[2]});
  }

  TrackerHit to_lcio2(EVENT::TrackerHit* h)
  {
    return TrackerHit(
      h->getCellID0(),
      h->getCellID1(),
      h->getTime(),
      h->getEDep(),
      h->getEDepError(),
      {h->getPosition()[0],h->getPosition()[1],h->getPosition()[2]});
  }

  ReconstructedParticle to_lcio2(EVENT::ReconstructedParticle* p)
  {
    return ReconstructedParticle(
      p->getType(),
      p->getEnergy(),
      {p->getMomentum()[0],p->getMomentum()[1],p->getMomentum()[2]},
      p->getCharge(),
      p->getMass());
  }

  Vertex to_lcio2(EVENT::Vertex* p)
  {
    return Vertex(
      p->isPrimary(),
      p->getChi2(),
      p->getProbability(),
      {p->getPosition()[0],p->getPosition()[1],p->getPosition()[2]});
  }

  Track to_lcio2(EVENT::Track* t)
  {
    return Track(
      t->getChi2(),
      t->getNdf(),
      t->getdEdx(),
      t->getdEdxError(),
      t->getRadiusOfInnermostHit());
  }
}

#endif
