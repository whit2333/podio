#ifndef LCIO2TOOLS_HH
#define LCIO2TOOLS_HH

#include "EVENT/MCParticle.h"
#include "EVENT/Cluster.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/SimCalorimeterHit.h"
#include "EVENT/RawCalorimeterHit.h"
#include "EVENT/TrackerHit.h"
#include "EVENT/TrackerRawData.h"
#include "EVENT/SimTrackerHit.h"
//#include "EVENT/RawTrackerHit.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/Vertex.h"

#include "IMPL/MCParticleImpl.h"
#include "IMPL/ClusterImpl.h"
#include "IMPL/CalorimeterHitImpl.h"
#include "IMPL/SimCalorimeterHitImpl.h"
#include "IMPL/RawCalorimeterHitImpl.h"
#include "IMPL/TrackerHitImpl.h"
#include "IMPL/TrackerRawDataImpl.h"
#include "IMPL/SimTrackerHitImpl.h"
//#include "IMPL/RawTrackerHitImpl.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/VertexImpl.h"
#include "IMPL/TrackImpl.h"

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
#include "lcio2/TrackerRawData.h"
#include "lcio2/Track.h"


namespace lcio2 {

  /** @defgroup  legacyhelpers Legacy helpers
   *  Helper functions for converting between legacy-lcio and lcio2
   */

  /** @defgroup tolcio2  Conversions to lcio2
   *  @ingroup legacyhelpers
   *  Conversions to lcio2 from legacy-lcio.
   *  @{
   */

  MCParticle            to_lcio2(EVENT::MCParticle* mcp);
  Cluster               to_lcio2(EVENT::Cluster* c);
  SimCalorimeterHit     to_lcio2(EVENT::SimCalorimeterHit* c);
  RawCalorimeterHit     to_lcio2(EVENT::RawCalorimeterHit* c);
  CalorimeterHit        to_lcio2(EVENT::CalorimeterHit* c);
  SimTrackerHit         to_lcio2(EVENT::SimTrackerHit* h);
  TrackerHit            to_lcio2(EVENT::TrackerHit* h);
  ReconstructedParticle to_lcio2(EVENT::ReconstructedParticle* p);
  Vertex                to_lcio2(EVENT::Vertex* p);
  Track                 to_lcio2(EVENT::Track* t);
  TrackerRawData        to_lcio2(EVENT::TrackerRawData* t);

  /** @} */

  /** @defgroup fromlcio2 Conversion to legacy-lcio
   *  @ingroup legacyhelpers
   *  Conversions from lcio2 to legacy-lcio.
   *  TODO: Finish legacy conversions. 
   *  @{
   */
  IMPL::MCParticleImpl*            from_lcio2(const MCParticle&            p);
  IMPL::ClusterImpl*               from_lcio2(const Cluster&               c);
  IMPL::SimCalorimeterHitImpl*     from_lcio2(const SimCalorimeterHit&     c);
  IMPL::RawCalorimeterHitImpl*     from_lcio2(const RawCalorimeterHit&     c);
  IMPL::CalorimeterHitImpl*        from_lcio2(const CalorimeterHit&        c);
  IMPL::SimTrackerHitImpl*         from_lcio2(const SimTrackerHit&         h);
  IMPL::TrackerHitImpl*            from_lcio2(const TrackerHit&            h);
  IMPL::ReconstructedParticleImpl* from_lcio2(const ReconstructedParticle& p);
  IMPL::VertexImpl*                from_lcio2(const Vertex&                p);
  IMPL::TrackImpl*                 from_lcio2(const Track&                 t);
  IMPL::TrackerRawDataImpl*        from_lcio2(const TrackerRawData&        t);
  /** @} */

}

#endif
