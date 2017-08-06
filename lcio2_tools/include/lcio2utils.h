#ifndef LCIO2TOOLS_HH
#define LCIO2TOOLS_HH

#include "EVENT/MCParticle.h"

#include "lcio2/MCParticleCollection.h"
#include "lcio2/MCParticle.h"

namespace lcio2 {

  MCParticle to_lcio2(EVENT::MCParticle* mcp) {
    return MCParticle(
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



}

#endif
