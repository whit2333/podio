#include "lcio2utils.h"
#include "IMPL/MCParticleImpl.h"

namespace lcio2 {

  MCParticle to_lcio2(EVENT::MCParticle* mcp)
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
      {0,0,0},
      c->getShape(),
      {0}, // FIXME: what are the weights used for?
      c->getSubdetectorEnergies()
      );
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

  TrackerRawData to_lcio2(EVENT::TrackerRawData* t)
  {
  TrackerRawData(
      t->getCellID0(),
      t->getCellID1(),
      0,//t->getChannelID(), // channel not used legacy-LCIO?
      t->getTime(),
      0);// adc = 0 ///FIXME
  }
  //______________________________________________________________________

  IMPL::MCParticleImpl*            from_lcio2(const MCParticle& h)
  {
    IMPL::MCParticleImpl* p = new IMPL::MCParticleImpl();
    p->setPDG(h.pdg());
    p->setGeneratorStatus(h.genstatus());
    p->setVertex(h.vertex().data());
    p->setCharge(h.charge());
    p->setMass(h.mass());
    p->setTime(h.time());
    if(h.endpointSet()){
      p->setEndpoint(h.endpoint().data());
    }
    p->setMomentum(h.momentum().data());
    return p;
  }

  IMPL::ClusterImpl*               from_lcio2(const Cluster&               c)
  {
    IMPL::ClusterImpl* p = new IMPL::ClusterImpl();
    p->setEnergy( c.energy() );
    p->setPosition( c.position().data() );
    p->setPositionError( c.positionError().data() );
    p->setITheta(  c.theta() );
    p->setIPhi(  c.phi() );
    p->setDirectionError( c.directionError().data() );
    p->setShape(c.shape());
    //p->setWeight(c.weight());
    p->subdetectorEnergies() = c.subdetectorEnergies();
    return p;
  }

  IMPL::SimCalorimeterHitImpl*     from_lcio2(const SimCalorimeterHit&     h)
  {
    IMPL::SimCalorimeterHitImpl* p = new IMPL::SimCalorimeterHitImpl();
    p->setCellID0(h.cellID0());
    p->setCellID1(h.cellID1());
    p->setEnergy(h.energy());
    //p->setTimeCont(h.time());
    p->setPosition(h.position().data());
    return p;
  }

  IMPL::RawCalorimeterHitImpl*     from_lcio2(const RawCalorimeterHit&     c)
  {
    IMPL::RawCalorimeterHitImpl* p = new IMPL::RawCalorimeterHitImpl();
    p->setCellID0(c.cellID0());
    p->setCellID1(c.cellID1());
    p->setAmplitude(c.amplitude());
    p->setTimeStamp(c.timeStamp());
    return p;
  }

  IMPL::CalorimeterHitImpl*        from_lcio2(const CalorimeterHit&        c)
  {
    IMPL::CalorimeterHitImpl* p = new IMPL::CalorimeterHitImpl();
    p->setCellID0(c.cellID0());
    p->setCellID1(c.cellID1());
    p->setEnergy(c.energy());
    p->setTime(c.time());
    p->setPosition(c.position().data());
    p->setType(c.type());
    return p;
  }

  IMPL::SimTrackerHitImpl*         from_lcio2(const SimTrackerHit&         h)
  {
    IMPL::SimTrackerHitImpl* p = new IMPL::SimTrackerHitImpl();
    p->setCellID0(h.cellID0());
    p->setCellID1(h.cellID1());
    p->setPosition( h.position().data());
    p->setTime(h.time());
    p->setPathLength(h.pathLength());
    p->setEDep(h.EDep());
    p->setMomentum(h._p().data());
    return p;
  }

  IMPL::TrackerHitImpl*            from_lcio2(const TrackerHit& h)
  {
    IMPL::TrackerHitImpl* p = new IMPL::TrackerHitImpl();
    p->setCellID0(h.cellID0());
    p->setCellID1(h.cellID1());
    p->setTime(h.time());
    p->setEDep(h.EDep());
    p->setEDepError(h.EDepError());
    p->setPosition( h.position().data());
    return p;
  }

  IMPL::ReconstructedParticleImpl* from_lcio2(const ReconstructedParticle& h)
  {
    IMPL::ReconstructedParticleImpl* p = new IMPL::ReconstructedParticleImpl();
    p->setType(h.type());
    p->setEnergy(h.energy());
    p->setMomentum(h.momentum().data());
    p->setCharge(h.charge());
    p->setMass(h.mass());
    return p;
  }

  IMPL::VertexImpl*                from_lcio2(const Vertex& h)
  {
    IMPL::VertexImpl* p = new IMPL::VertexImpl();
      p->setPrimary(h.primary());
      p->setChi2(h.chi2());
      p->setProbability(h.probability());
      p->setPosition(h.position().data());
    return p;
  }

  IMPL::TrackImpl*                 from_lcio2(const Track& h)
  {
    IMPL::TrackImpl* p = new IMPL::TrackImpl();
    p->setChi2(h.chi2());
    p->setNdf(h.ndf());
    p->setdEdx(h.dEdx());
    p->setdEdxError(h.dEdxError());
    p->setRadiusOfInnermostHit(h.radiusOfInnermostHit());
    return p;
  }

  IMPL::TrackerRawDataImpl*        from_lcio2(const TrackerRawData& h)
  {
    IMPL::TrackerRawDataImpl* p = new IMPL::TrackerRawDataImpl();
    p->setCellID0(h.cellID0());
    p->setCellID1(h.cellID1());
    p->setTime(h.time());
    return p;
  }

}

