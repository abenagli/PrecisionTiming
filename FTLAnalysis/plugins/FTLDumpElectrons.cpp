#ifndef _FTL_DUMP_ELECTRONS_
#define _FTL_DUMP_ELECTRONS_

#include "TMath.h"

#include "FWCore/Utilities/interface/BranchType.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/Provenance.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"
#include "DataFormats/FTLRecHit/interface/FTLUncalibratedRecHit.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHit.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"
#include "DataFormats/FTLRecHit/interface/FTLClusterCollections.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementError.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "RecoMTD/TransientTrackingRecHit/interface/MTDTransientTrackingRecHitBuilder.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderWithPropagator.h"
#include "RecoTracker/TransientTrackingRecHit/interface/Traj2TrackHits.h"
#include "TrackingTools/TrackRefitter/interface/TrackTransformer.h"

#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "Geometry/CommonTopologies/interface/Topology.h"
#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeomDetUnit.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"
#include "DataFormats/GeometrySurface/interface/BoundSurface.h"
#include "DataFormats/GeometrySurface/interface/MediumProperties.h"
#include "DataFormats/GeometrySurface/interface/TrapezoidalPlaneBounds.h"

#include "RecoMTD/DetLayers/interface/MTDDetLayerGeometry.h"
#include "RecoMTD/DetLayers/interface/MTDTrayBarrelLayer.h"
#include "RecoMTD/DetLayers/interface/MTDDetTray.h"
#include "RecoMTD/DetLayers/interface/MTDRingForwardDoubleLayer.h"
#include "RecoMTD/DetLayers/interface/MTDDetRing.h"
#include "RecoMTD/Records/interface/MTDRecoGeometryRecord.h"

#include "PrecisionTiming/FTLAnalysis/interface/FTLElectronsTree.h"

#include <memory>
using namespace std;



class FTLDumpElectrons : public edm::EDAnalyzer
{
public:
  explicit FTLDumpElectrons(const edm::ParameterSet& pSet);
  ~FTLDumpElectrons() {};
  
  //---utils
  
  //---methods
  virtual void beginJob() override {};
  virtual void analyze(edm::Event const&, edm::EventSetup const&) override;
  virtual void endJob() override {};
  
  std::string PrintPosition(const GlobalPoint& gp);
  std::string PrintPosition(const LocalPoint& lp);
  
private:
  const MTDGeometry* mtdGeometry_;
  
  //---inputs
  edm::Handle<reco::GenParticleCollection> genParticlesHandle_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::Handle<vector<SimVertex> > genVtxHandle_;
  edm::EDGetTokenT<vector<SimVertex> > genVtxToken_;
  
  edm::Handle<std::vector<PSimHit> > simHitsBTLHandle_;
  edm::EDGetTokenT<std::vector<PSimHit> > simHitsBTLToken_;
  edm::Handle<FTLRecHitCollection> recHitsBTLHandle_;
  edm::EDGetTokenT<FTLRecHitCollection> recHitsBTLToken_;
  edm::Handle<FTLClusterCollection> clustersBTLHandle_;
  edm::EDGetTokenT<FTLClusterCollection> clustersBTLToken_;    
  
  edm::Handle<std::vector<PSimHit> > simHitsETLHandle_;
  edm::EDGetTokenT<std::vector<PSimHit> > simHitsETLToken_;
  edm::Handle<FTLRecHitCollection> recHitsETLHandle_;
  edm::EDGetTokenT<FTLRecHitCollection> recHitsETLToken_;    
  edm::Handle<FTLClusterCollection> clustersETLHandle_;
  edm::EDGetTokenT<FTLClusterCollection> clustersETLToken_;    
  
  edm::EDGetTokenT<reco::GsfElectronCollection> barrelElectronsToken_;
  edm::Handle<reco::GsfElectronCollection> barrelElectronsHandle_;
  edm::EDGetTokenT<reco::GsfElectronCollection> endcapElectronsToken_;
  edm::Handle<reco::GsfElectronCollection> endcapElectronsHandle_;
  
  edm::EDGetTokenT<edm::ValueMap<float> > barrelElectronMVAToken_;
  edm::Handle<edm::ValueMap<float> > barrelElectronMVAHandle_;
  edm::EDGetTokenT<edm::ValueMap<float> > endcapElectronMVAToken_;
  edm::Handle<edm::ValueMap<float> > endcapElectronMVAHandle_;
  
  edm::ESHandle<TransientTrackBuilder> ttrackBuilder;
  
  //---options
  BTLDetId::CrysLayout crysLayout_;
  double track_hit_DRMax_;
  double track_hit_distMax_;
  bool verbosity_;
  
  //---outputs
  FTLElectronsTree outTree_;
  edm::Service<TFileService> fs_;  
};



FTLDumpElectrons::FTLDumpElectrons(const edm::ParameterSet& pSet):
  genParticlesToken_(consumes<reco::GenParticleCollection>(pSet.getUntrackedParameter<edm::InputTag>("genParticlesTag"))),
  genVtxToken_(consumes<vector<SimVertex> >(pSet.getUntrackedParameter<edm::InputTag>("genVtxTag"))),
  simHitsBTLToken_(consumes<std::vector<PSimHit> >(pSet.getUntrackedParameter<edm::InputTag>("simHitsBTLTag"))),
  recHitsBTLToken_(consumes<FTLRecHitCollection>(pSet.getUntrackedParameter<edm::InputTag>("recHitsBTLTag"))),
  clustersBTLToken_(consumes<FTLClusterCollection>(pSet.getUntrackedParameter<edm::InputTag>("clustersBTLTag"))),
  simHitsETLToken_(consumes<std::vector<PSimHit> >(pSet.getUntrackedParameter<edm::InputTag>("simHitsETLTag"))),
  recHitsETLToken_(consumes<FTLRecHitCollection>(pSet.getUntrackedParameter<edm::InputTag>("recHitsETLTag"))),
  clustersETLToken_(consumes<FTLClusterCollection>(pSet.getUntrackedParameter<edm::InputTag>("clustersETLTag"))),
  barrelElectronsToken_(consumes<reco::GsfElectronCollection>(pSet.getUntrackedParameter<edm::InputTag>("barrelElectronsTag"))),
  endcapElectronsToken_(consumes<reco::GsfElectronCollection>(pSet.getUntrackedParameter<edm::InputTag>("endcapElectronsTag"))),
  barrelElectronMVAToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("barrelElectronMVATag"))),
  endcapElectronMVAToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("endcapElectronMVATag"))),
  crysLayout_((BTLDetId::CrysLayout)(pSet.getUntrackedParameter<int>("crysLayout"))),
  track_hit_DRMax_(pSet.getParameter<double>("track_hit_DRMax")),
  track_hit_distMax_(pSet.getParameter<double>("track_hit_distMax")),
  verbosity_(pSet.getParameter<bool>("verbosity"))
{
  outTree_ = FTLElectronsTree(pSet.getUntrackedParameter<string>("treeName").c_str(), "FTLHits tree for FTL studies");
}



void FTLDumpElectrons::analyze(edm::Event const& event, edm::EventSetup const& setup)
{
  outTree_.Reset();
  
  
  //---get the MTD geometry
  edm::ESHandle<MTDGeometry> geoHandle;
  setup.get<MTDDigiGeometryRecord>().get(geoHandle);
  mtdGeometry_ = geoHandle.product();
  
  edm::ESHandle<MTDDetLayerGeometry> layerGeo;
  setup.get<MTDRecoGeometryRecord>().get(layerGeo);
  
  
  //--- get the B field
  edm::ESHandle<MagneticField> theField;
  setup.get<IdealMagneticFieldRecord>().get(theField);
  
  
  //--- load gen particles
  event.getByToken(genParticlesToken_, genParticlesHandle_);
  auto genParticles = *genParticlesHandle_.product();
  
  event.getByToken(genVtxToken_, genVtxHandle_);    
  const SimVertex* genPV = NULL;
  if(genVtxHandle_.isValid()) genPV = &(genVtxHandle_.product()->at(0));
  auto simVertices = *genVtxHandle_.product();
  
  //---load BTL hits
  event.getByToken(simHitsBTLToken_, simHitsBTLHandle_);
  auto simHitsBTL = *simHitsBTLHandle_.product();
  
  event.getByToken(recHitsBTLToken_, recHitsBTLHandle_);
  auto recHitsBTL = FTLRecHitCollection();
  if(recHitsBTLHandle_.isValid()) recHitsBTL = *recHitsBTLHandle_.product();
  
  event.getByToken(clustersBTLToken_, clustersBTLHandle_);
  auto clustersBTL = *clustersBTLHandle_.product();
  
  
  //---load ETL hits
  event.getByToken(simHitsETLToken_, simHitsETLHandle_);
  auto simHitsETL = *simHitsETLHandle_.product();
  
  //---load the FTL collection if present in the EventContent (avoid crash with standard geometry)
  event.getByToken(recHitsETLToken_, recHitsETLHandle_);
  auto recHitsETL = FTLRecHitCollection();
  if(recHitsETLHandle_.isValid()) recHitsETL = *recHitsETLHandle_.product();
  
  event.getByToken(clustersETLToken_, clustersETLHandle_);
  auto clustersETL = *clustersETLHandle_.product();
  
  
  //---load electrons
  event.getByToken(barrelElectronsToken_,barrelElectronsHandle_);
  auto barrelElectrons = *barrelElectronsHandle_.product();
  
  event.getByToken(endcapElectronsToken_,endcapElectronsHandle_);
  auto endcapElectrons = *endcapElectronsHandle_.product();
  
  event.getByToken(barrelElectronMVAToken_,barrelElectronMVAHandle_);
  auto barrelElectronMVA = *barrelElectronMVAHandle_.product();
  
  event.getByToken(endcapElectronMVAToken_,endcapElectronMVAHandle_);
  auto endcapElectronMVA = *endcapElectronMVAHandle_.product();
  
  
  //---load tracks
  setup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttrackBuilder);
  
  
  
  
  
  
  //---fill the tree - gen vertex
  if( genPV )
  {
    outTree_.genVtx_x = genPV->position().x();
    outTree_.genVtx_y = genPV->position().y();
    outTree_.genVtx_z = genPV->position().z();
    outTree_.genVtx_t = genPV->position().t()*1E9; //ns                                                                                                                                                                                    
  }
  
  
  
  //--- fill the tree - barrel electrons
  int idx=0;
  for(edm::View<reco::GsfElectron>::size_type iEle = 0; iEle < barrelElectrons.size(); ++iEle)
  {
    auto electron = barrelElectrons.at(iEle);
    edm::Ref<reco::GsfElectronCollection> electronRef(barrelElectronsHandle_,iEle);
    auto trackRef = electron.gsfTrack();
    if( trackRef.isNull() ) continue;
    
    
    // skip neutrals
    if( electron.charge() == 0 ) continue;
    if( electron.pt() < 0.7 ) continue;
    
    // match with gen particles
    float DRMin = 999999;
    int genPdgId = 0;
    float genEta = -999.;
    float genPhi = -999.;
    float genPt = -999.;
    for(auto& genPart : genParticles)
    {
      if( genPart.status() != 1 ) continue;
      if( genPart.charge() == 0 ) continue;
      
      float Deta = electron.eta()-genPart.eta();
      float Dphi = deltaPhi(electron.phi(),genPart.phi());
      float DR   = sqrt(Deta*Deta+Dphi*Dphi);
      
      if( DR < DRMin )
      {
        DRMin = DR;
        
        genPdgId = genPart.pdgId();
        genEta   = genPart.eta();
        genPhi   = genPart.phi();
        genPt    = genPart.pt();
      }
    }
    
    if( genPV )
    {
      outTree_.electrons_mcMatch_genVtx_x -> push_back(genPV->position().x());
      outTree_.electrons_mcMatch_genVtx_y -> push_back(genPV->position().y());
      outTree_.electrons_mcMatch_genVtx_z -> push_back(genPV->position().z());
      outTree_.electrons_mcMatch_genVtx_t -> push_back(genPV->position().t()*1E9); //ns
    }
    else
    {
      outTree_.electrons_mcMatch_genVtx_x -> push_back(-999.);
      outTree_.electrons_mcMatch_genVtx_y -> push_back(-999.);
      outTree_.electrons_mcMatch_genVtx_z -> push_back(-999.);
      outTree_.electrons_mcMatch_genVtx_t -> push_back(-999.);
    }
    
    outTree_.electrons_idx -> push_back(idx);
    outTree_.electrons_pt -> push_back(electron.pt());
    outTree_.electrons_eta -> push_back(electron.eta());
    outTree_.electrons_phi -> push_back(electron.phi());
    outTree_.electrons_x -> push_back(electron.vx());
    outTree_.electrons_y -> push_back(electron.vy());
    outTree_.electrons_z -> push_back(electron.vz());
    outTree_.electrons_t -> push_back(trackRef->t0());
    outTree_.electrons_energy -> push_back(electron.energy());
    outTree_.electrons_hasMTD -> push_back(trackRef->isTimeOk());
    outTree_.electrons_mcMatch_genPdgId -> push_back(genPdgId);
    outTree_.electrons_mcMatch_genPt -> push_back(genPt);
    outTree_.electrons_mcMatch_genEta -> push_back(genEta);
    outTree_.electrons_mcMatch_genPhi -> push_back(genPhi);
    outTree_.electrons_mcMatch_DR -> push_back(DRMin);
    outTree_.electrons_mva -> push_back(barrelElectronMVA[electronRef]);
    
    outTree_.electrons_isEB -> push_back(int(electron.isEB()));
    outTree_.electrons_isEE -> push_back(int(electron.isEE()));
    outTree_.electrons_fbrem -> push_back(electron.fbrem());
    outTree_.electrons_dEtaIn -> push_back(electron.deltaEtaSuperClusterTrackAtVtx());
    outTree_.electrons_dPhiIn -> push_back(electron.deltaPhiSuperClusterTrackAtVtx());
    outTree_.electrons_eop -> push_back(electron.eSuperClusterOverP());
    outTree_.electrons_sigmaEtaEta -> push_back(electron.sigmaEtaEta());
    outTree_.electrons_sigmaIetaIeta -> push_back(electron.full5x5_sigmaIetaIeta());
    outTree_.electrons_sigmaIphiIphi -> push_back(electron.full5x5_sigmaIphiIphi());
    outTree_.electrons_ecalEnergy -> push_back(electron.ecalEnergy());
    outTree_.electrons_e1x5 -> push_back(electron.full5x5_e1x5());
    outTree_.electrons_e2x5Max -> push_back(electron.full5x5_e2x5Max());
    outTree_.electrons_e5x5 -> push_back(electron.full5x5_e5x5());
    outTree_.electrons_r9 -> push_back(electron.full5x5_r9());
    outTree_.electrons_hcalOverEcal -> push_back(electron.hcalOverEcal());
    
    if( verbosity_ ) std::cout << "*** electron " << iEle << " / " << barrelElectrons.size() << "   pt: " << electron.pt() << "   eta: " << electron.eta() << "   phi: " << electron.phi() << std::endl;
    if( verbosity_ ) std::cout << "*** match with gen particle   DR: " << DRMin << "   gen pdgId: " << genPdgId << "   gen eta: " << genEta << "   gen phi: " << genPhi << "   genPt: " << genPt << std::endl;
    
    outTree_.matchedSimHits_n->resize(idx+1);
    outTree_.matchedRecHits_n->resize(idx+1);
    outTree_.matchedClusters_n->resize(idx+1);
    outTree_.matchedSimHits_idx->resize(idx+1);
    outTree_.matchedSimHits_det->resize(idx+1);
    outTree_.matchedSimHits_energy->resize(idx+1);
    outTree_.matchedSimHits_energyCorr->resize(idx+1);
    outTree_.matchedSimHits_time->resize(idx+1);
    outTree_.matchedSimHits_rr->resize(idx+1);
    outTree_.matchedSimHits_module->resize(idx+1);
    outTree_.matchedSimHits_modType->resize(idx+1);
    outTree_.matchedSimHits_crystal->resize(idx+1);
    outTree_.matchedSimHits_ieta->resize(idx+1);
    outTree_.matchedSimHits_iphi->resize(idx+1);
    outTree_.matchedSimHits_entry_local_x->resize(idx+1);
    outTree_.matchedSimHits_entry_local_y->resize(idx+1);
    outTree_.matchedSimHits_entry_local_z->resize(idx+1);
    outTree_.matchedSimHits_entry_global_R->resize(idx+1);
    outTree_.matchedSimHits_exit_local_x->resize(idx+1);
    outTree_.matchedSimHits_exit_local_y->resize(idx+1);
    outTree_.matchedSimHits_exit_local_z->resize(idx+1);
    outTree_.matchedSimHits_exit_global_R->resize(idx+1);
    outTree_.matchedSimHits_electron_Deta->resize(idx+1);
    outTree_.matchedSimHits_electron_Dphi->resize(idx+1);
    outTree_.matchedSimHits_electron_DR->resize(idx+1);
    outTree_.matchedSimHits_electron_Dz->resize(idx+1);
    outTree_.matchedSimHits_electron_RDphi->resize(idx+1);
    outTree_.matchedSimHits_electron_dist->resize(idx+1);
    outTree_.matchedRecHits_idx->resize(idx+1);
    outTree_.matchedRecHits_det->resize(idx+1);
    outTree_.matchedRecHits_energy->resize(idx+1);
    outTree_.matchedRecHits_energyCorr->resize(idx+1);
    outTree_.matchedRecHits_time->resize(idx+1);
    outTree_.matchedRecHits_rr->resize(idx+1);
    outTree_.matchedRecHits_module->resize(idx+1);
    outTree_.matchedRecHits_modType->resize(idx+1);
    outTree_.matchedRecHits_crystal->resize(idx+1);
    outTree_.matchedRecHits_ieta->resize(idx+1);
    outTree_.matchedRecHits_iphi->resize(idx+1);
    outTree_.matchedRecHits_local_x->resize(idx+1);
    outTree_.matchedRecHits_local_y->resize(idx+1);
    outTree_.matchedRecHits_local_z->resize(idx+1);
    outTree_.matchedRecHits_global_R->resize(idx+1);
    outTree_.matchedRecHits_electron_Deta->resize(idx+1);
    outTree_.matchedRecHits_electron_Dphi->resize(idx+1);
    outTree_.matchedRecHits_electron_DR->resize(idx+1);
    outTree_.matchedRecHits_electron_Dz->resize(idx+1);
    outTree_.matchedRecHits_electron_RDphi->resize(idx+1);
    outTree_.matchedRecHits_electron_dist->resize(idx+1);
    outTree_.matchedRecHits_sietaieta->resize(idx+1);
    outTree_.matchedRecHits_siphiiphi->resize(idx+1);
    outTree_.matchedClusters_idx->resize(idx+1);
    outTree_.matchedClusters_det->resize(idx+1);
    outTree_.matchedClusters_energy->resize(idx+1);
    outTree_.matchedClusters_energyCorr->resize(idx+1);
    outTree_.matchedClusters_time->resize(idx+1);
    outTree_.matchedClusters_rr->resize(idx+1);
    outTree_.matchedClusters_module->resize(idx+1);
    outTree_.matchedClusters_modType->resize(idx+1);
    outTree_.matchedClusters_crystal->resize(idx+1);
    outTree_.matchedClusters_ieta->resize(idx+1);
    outTree_.matchedClusters_iphi->resize(idx+1);
    outTree_.matchedClusters_size->resize(idx+1);
    outTree_.matchedClusters_size_x->resize(idx+1);
    outTree_.matchedClusters_size_y->resize(idx+1);
    outTree_.matchedClusters_local_x->resize(idx+1);
    outTree_.matchedClusters_local_y->resize(idx+1);
    outTree_.matchedClusters_local_z->resize(idx+1);
    outTree_.matchedClusters_global_R->resize(idx+1);
    outTree_.matchedClusters_electron_Deta->resize(idx+1);
    outTree_.matchedClusters_electron_Dphi->resize(idx+1);
    outTree_.matchedClusters_electron_DR->resize(idx+1);
    outTree_.matchedClusters_electron_Dz->resize(idx+1);
    outTree_.matchedClusters_electron_RDphi->resize(idx+1);
    outTree_.matchedClusters_electron_dist->resize(idx+1);
        
    //---get compatible layers/Dets
    std::vector<GlobalPoint> gp_ext;
    std::vector<LocalPoint> lp_ext;
    auto tTrack = ttrackBuilder->build(trackRef);
    TrajectoryStateOnSurface tsos = tTrack.outermostMeasurementState();
    float theMaxChi2 = 25.;
    float theNSigma = 5.;
    std::unique_ptr<MeasurementEstimator> theEstimator = std::make_unique<Chi2MeasurementEstimator>(theMaxChi2,theNSigma);
    SteppingHelixPropagator prop(theField.product(),anyDirection);
    
    //try BTL
    bool inBTL = false;
    const vector<const DetLayer*>& layersBTL = layerGeo->allBTLLayers();
    for(const DetLayer* ilay : layersBTL) 
    {
      pair<bool, TrajectoryStateOnSurface> comp = ilay->compatible(tsos,prop,*theEstimator);	
      if (!comp.first) continue;
      vector<DetLayer::DetWithState> compDets = ilay->compatibleDets(tsos,prop,*theEstimator);
      for( const auto& detWithState : compDets ) 
      {
        const auto& det = detWithState.first;
        const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
        const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());
        
        gp_ext.push_back(detWithState.second.globalPosition());
        lp_ext.push_back(topo.moduleToPixelLocalPoint(det->toLocal(gp_ext.back())));
        if( !inBTL )
        {
          outTree_.electrons_eta_atBTL -> push_back(gp_ext.back().eta());
          outTree_.electrons_phi_atBTL -> push_back(gp_ext.back().phi());
          outTree_.electrons_local_x_atBTL -> push_back(lp_ext.back().x());
          outTree_.electrons_local_y_atBTL -> push_back(lp_ext.back().y());
          outTree_.electrons_local_z_atBTL -> push_back(lp_ext.back().z());
          outTree_.electrons_global_R_atBTL -> push_back(sqrt(gp_ext.back().perp2()));
          inBTL = true;
          
          if( verbosity_ )
            std::cout << ">>> extrapolated point in BTL:   global: " << gp_ext.back() << "   local: " << lp_ext.back() << std::endl;
        }
      }
    }
    
    //try ETL
    bool inETL = false;
    const vector<const DetLayer*>& layersETL = layerGeo->allETLLayers();
    for(const DetLayer* ilay : layersETL) 
    {
      const BoundDisk& disk = static_cast<const MTDRingForwardDoubleLayer*>(ilay)->specificSurface();
      const double diskZ = disk.position().z();
      if( tsos.globalPosition().z() * diskZ < 0 ) continue; // only propagate to the disk that's on the same side
      pair<bool, TrajectoryStateOnSurface> comp = ilay->compatible(tsos,prop,*theEstimator);	
      if (!comp.first) continue;
      vector<DetLayer::DetWithState> compDets = ilay->compatibleDets(tsos,prop,*theEstimator);
      for( const auto& detWithState : compDets ) 
      {
        const auto& det = detWithState.first;
        const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
        const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());
        
        gp_ext.push_back(detWithState.second.globalPosition());
        lp_ext.push_back(topo.moduleToPixelLocalPoint(detWithState.first->toLocal(gp_ext.back())));
        if( !inETL )
        {
          outTree_.electrons_eta_atETL -> push_back(gp_ext.back().eta());
          outTree_.electrons_phi_atETL -> push_back(gp_ext.back().phi());
          outTree_.electrons_local_x_atETL -> push_back(lp_ext.back().x());
          outTree_.electrons_local_y_atETL -> push_back(lp_ext.back().y());
          outTree_.electrons_local_z_atETL -> push_back(lp_ext.back().z());
          outTree_.electrons_global_R_atETL -> push_back(sqrt(gp_ext.back().perp2()));
          inETL=true;
          
          if( verbosity_ )
            std::cout << ">>> extrapolated point in ETL:   global: " << gp_ext.back() << "   local: " << lp_ext.back() << std::endl;
        }
      }
    }
    
    if( !inBTL )
    {
      outTree_.electrons_eta_atBTL -> push_back(-999.);
      outTree_.electrons_phi_atBTL -> push_back(-999.);
      outTree_.electrons_local_x_atBTL -> push_back(-999.);
      outTree_.electrons_local_y_atBTL -> push_back(-999.);
      outTree_.electrons_local_z_atBTL -> push_back(-999.);
      outTree_.electrons_global_R_atBTL -> push_back(-999.);
    }
    
    if( !inETL )
    {
      outTree_.electrons_eta_atETL -> push_back(-999.);
      outTree_.electrons_phi_atETL -> push_back(-999.);
      outTree_.electrons_local_x_atETL -> push_back(-999.);
      outTree_.electrons_local_y_atETL -> push_back(-999.);
      outTree_.electrons_local_z_atETL -> push_back(-999.);
      outTree_.electrons_global_R_atETL -> push_back(-999.);
    }
    if( verbosity_ ) std::cout << "---" << std::endl;

    
    
    //---get associated BTL simHits
    if( verbosity_ ) std::cout << "*** simHits - n tot: " << simHitsBTL.size() << std::endl;
    int simHitIt = 0;
    for(auto simHit : simHitsBTL)
    {
      BTLDetId id = simHit.detUnitId();
      DetId geoId = id.geographicalId( crysLayout_ );
      const auto& det = mtdGeometry_ -> idToDet(geoId);
      const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
      const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());
      
      double energy = simHit.energyLoss()*1000.;
      double time   = simHit.tof();
      
      if ((time)<0 || (time)>25) continue;
      int RR = id.mtdRR();
      int module = id.module();
      int modType = id.modType();
      int crystal = id.crystal();
      int ieta = id.ieta(crysLayout_);
      int iphi = id.iphi(crysLayout_);
      
      LocalPoint lp_entry(   simHit.entryPoint().x()/10.,   simHit.entryPoint().y()/10.,   simHit.entryPoint().z()/10.);
      LocalPoint lp_mid  (simHit.localPosition().x()/10.,simHit.localPosition().y()/10.,simHit.localPosition().z()/10.);
      LocalPoint lp_exit (    simHit.exitPoint().x()/10.,    simHit.exitPoint().y()/10.,    simHit.exitPoint().z()/10.);
      GlobalPoint gp_entry = det->toGlobal(topo.pixelToModuleLocalPoint(lp_entry,id.row(topo.nrows()),id.column(topo.nrows())));
      GlobalPoint gp_mid   = det->toGlobal(topo.pixelToModuleLocalPoint(lp_mid,id.row(topo.nrows()),id.column(topo.nrows())));
      GlobalPoint gp_exit  = det->toGlobal(topo.pixelToModuleLocalPoint(lp_exit,id.row(topo.nrows()),id.column(topo.nrows())));
      
      float eta = gp_mid.eta();
      float phi = gp_mid.phi();
      
      int closestPoint=-1;
      float minDist=999;
      for (unsigned int ic=0; ic<gp_ext.size();++ic)
      {
        GlobalPoint diff(gp_ext[ic].x()-gp_mid.x(),gp_ext[ic].y()-gp_mid.y(),gp_ext[ic].z()-gp_mid.z());
        if (diff.mag()<minDist)
        {
          closestPoint=ic;
          minDist=diff.mag();
        }
      }
      
      if (closestPoint == -1)
        continue;
      
      GlobalPoint gp_track = gp_ext[closestPoint];      
      
      float Deta  = gp_track.mag() > 0. ? eta-gp_track.eta()           : -999.;
      float Dphi  = gp_track.mag() > 0. ? deltaPhi(phi,gp_track.phi()) : -999.;
      float DR    = gp_track.mag() > 0. ? sqrt(Deta*Deta+Dphi*Dphi)    : -999.;
      float Dz    = gp_track.mag() > 0. ? gp_track.z()-gp_mid.z()      : -999.;
      float RDphi = gp_track.mag() > 0. ? sqrt(gp_track.perp2())*Dphi  : -999.;
      float dist  = gp_track.mag() > 0. ? (gp_mid-gp_track).mag()      : -999.;
      if( DR < 0.05 && DR > 0. )
      {
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ":   energy: " << energy << " MeV   time: " << time << " ns"
                                   << "   RR: " << RR << "   module: " << module << "   modType: " << modType << "   crystal: " << crystal
                                   << "   ieta: " << ieta << "   iphi: " << iphi << "   row: " << id.row(topo.nrows()) << "   column: " << id.column(topo.nrows()) << std::endl;
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ":   entryPoint:   local: " << PrintPosition(lp_entry)        << "   global: " << PrintPosition(gp_entry) << std::endl;
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ":     midPoint:   local: " << PrintPosition(lp_mid)          << "   global: " << PrintPosition(gp_mid)   << std::endl;
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ":    exitPoint:   local: " << PrintPosition(lp_exit)         << "   global: " << PrintPosition(gp_exit)  << std::endl;
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ":  detPosition:  global: " << PrintPosition(det->position()) << std::endl;
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ": hit position:   local: " << PrintPosition(lp_mid)          << "   global: " << PrintPosition(gp_mid) << "   DR: " << DR << "   dist: " << (gp_track-gp_mid).mag() << std::endl;
      }
      
      if( DR > track_hit_DRMax_ || dist > track_hit_distMax_ ) continue;
      if( gp_track.mag() <= 0. ) continue;
      
      outTree_.matchedSimHits_n->at(idx) += 1;
      
      outTree_.matchedSimHits_idx->at(idx).push_back(idx);
      outTree_.matchedSimHits_det->at(idx).push_back(1);
      outTree_.matchedSimHits_energy->at(idx).push_back(energy);
      outTree_.matchedSimHits_energyCorr->at(idx).push_back(energy*fabs(sin(trackRef->theta())));
      outTree_.matchedSimHits_time->at(idx).push_back(time);
      outTree_.matchedSimHits_rr->at(idx).push_back(RR);
      outTree_.matchedSimHits_module->at(idx).push_back(module);
      outTree_.matchedSimHits_modType->at(idx).push_back(modType);
      outTree_.matchedSimHits_crystal->at(idx).push_back(crystal);
      outTree_.matchedSimHits_ieta->at(idx).push_back(ieta);
      outTree_.matchedSimHits_iphi->at(idx).push_back(iphi);
      outTree_.matchedSimHits_entry_local_x->at(idx).push_back(lp_entry.x());
      outTree_.matchedSimHits_entry_local_y->at(idx).push_back(lp_entry.y());
      outTree_.matchedSimHits_entry_local_z->at(idx).push_back(lp_entry.z());
      outTree_.matchedSimHits_entry_global_R->at(idx).push_back(sqrt(gp_entry.perp2()));
      outTree_.matchedSimHits_exit_local_x->at(idx).push_back(lp_exit.x());
      outTree_.matchedSimHits_exit_local_y->at(idx).push_back(lp_exit.y());
      outTree_.matchedSimHits_exit_local_z->at(idx).push_back(lp_exit.z());
      outTree_.matchedSimHits_exit_global_R->at(idx).push_back(sqrt(gp_exit.perp2()));
      outTree_.matchedSimHits_electron_Deta->at(idx).push_back(fabs(Deta));
      outTree_.matchedSimHits_electron_Dphi->at(idx).push_back(fabs(Dphi));
      outTree_.matchedSimHits_electron_DR->at(idx).push_back(DR);
      outTree_.matchedSimHits_electron_Dz->at(idx).push_back(Dz);
      outTree_.matchedSimHits_electron_RDphi->at(idx).push_back(RDphi);
      outTree_.matchedSimHits_electron_dist->at(idx).push_back(dist);
      
      ++simHitIt;
    }
    if( verbosity_ ) std::cout << "---" << std::endl;
    
    
    //---find associated BTL recHits
    float sieie=0, sipip=0;
    float ss_hit_count=0;
    int recHitIt = 0;
    if( verbosity_ ) std::cout << "*** recHits - tot: " << recHitsBTL.size() << std::endl;
    for(auto recHit : recHitsBTL)
    {
      BTLDetId id = recHit.id();
      DetId geoId = id.geographicalId( crysLayout_ );
      const auto& det = mtdGeometry_ -> idToDet(geoId);
      const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
      const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());    
      
      double energy = recHit.energy();
      double time   = recHit.time();
      
      MeasurementPoint mp(recHit.row()+0.5f,recHit.column()+0.5f); //to get the center you need +0.5
      LocalPoint lp = topo.localPosition(mp);
      GlobalPoint gp = det->toGlobal(lp);

      float eta = gp.eta();
      float phi = gp.phi();
      
      int RR = id.mtdRR();
      int module = id.module();
      int modType = id.modType();
      int crystal = id.crystal();
      int ieta = id.ieta(crysLayout_);
      int iphi = id.iphi(crysLayout_);
      
      int closestPoint=-1;
      float minDist=999;
      for (unsigned int ic=0; ic<gp_ext.size();++ic)
      {
        GlobalPoint diff(gp_ext[ic].x()-gp.x(),gp_ext[ic].y()-gp.y(),gp_ext[ic].z()-gp.z());
        if( diff.mag() < minDist )
        {
          closestPoint=ic;
          minDist=diff.mag();
        }
      }
      
      if (closestPoint == -1) continue;
      
      GlobalPoint gp_track = gp_ext[closestPoint];      
      
      float Deta  = gp_track.mag() > 0. ? eta-gp_track.eta()           : -999.;
      float Dphi  = gp_track.mag() > 0. ? deltaPhi(phi,gp_track.phi()) : -999.;
      float DR    = gp_track.mag() > 0. ? sqrt(Deta*Deta+Dphi*Dphi)    : -999.;
      float Dz    = gp_track.mag() > 0. ? gp_track.z()-gp.z()          : -999.;
      float RDphi = gp_track.mag() > 0. ? sqrt(gp_track.perp2())*Dphi  : -999.;
      float dist  = gp_track.mag() > 0. ? (gp-gp_track).mag()          : -999.;
      
      // if( DR < 0.2 && DR > 0. )
      {
        if( verbosity_ ) std::cout << ">>> " << recHitIt << ":   energy: " << energy << " MeV   time: " << time << " ns"
                                   << "   RR: " << RR << "   module: " << module << "   modType: " << modType << "   crystal: " << crystal
                                   << "   ieta: " << ieta << "   iphi: " << iphi << "   row: " << recHit.row() << "- " << id.row(topo.nrows()) << "   column: " << recHit.column() << std::endl;
        if( verbosity_ ) std::cout << ">>> " << recHitIt << ":  detPosition:  global: " << PrintPosition(det->position()) << "   DR: " << DR << "   dist: " << (gp_track-det->position()).mag() << std::endl;
        if( verbosity_ ) std::cout << ">>> " << recHitIt << ": hit position:   local: " << PrintPosition(lp) << "   global: " << PrintPosition(gp) << "DR: " << DR << "   dist: " << (gp_track-gp).mag() << std::endl;
        if( verbosity_ ) std::cout << ">>> " << recHitIt << ": extrapolated track " << abs(recHit.row()-int()) << " : " << PrintPosition(gp_track) << std::endl;
      }
      
      if( DR > track_hit_DRMax_ || dist > track_hit_distMax_ ) continue;
      if( gp_track.mag() <= 0. ) continue;
      
      outTree_.matchedRecHits_n->at(idx) += 1;
      
      outTree_.matchedRecHits_idx->at(idx).push_back(idx);
      outTree_.matchedRecHits_det->at(idx).push_back(1);
      outTree_.matchedRecHits_energy->at(idx).push_back(energy);
      outTree_.matchedRecHits_energyCorr->at(idx).push_back(energy*fabs(sin(trackRef->theta())));
      outTree_.matchedRecHits_time->at(idx).push_back(time);
      outTree_.matchedRecHits_rr->at(idx).push_back(RR);
      outTree_.matchedRecHits_module->at(idx).push_back(module);
      outTree_.matchedRecHits_modType->at(idx).push_back(modType);
      outTree_.matchedRecHits_crystal->at(idx).push_back(crystal);
      outTree_.matchedRecHits_ieta->at(idx).push_back(ieta);
      outTree_.matchedRecHits_iphi->at(idx).push_back(iphi);
      outTree_.matchedRecHits_local_x->at(idx).push_back(lp.x());
      outTree_.matchedRecHits_local_y->at(idx).push_back(lp.y());
      outTree_.matchedRecHits_local_z->at(idx).push_back(lp.z());
      outTree_.matchedRecHits_global_R->at(idx).push_back(sqrt(gp.perp2()));
      outTree_.matchedRecHits_electron_Deta->at(idx).push_back(fabs(Deta));
      outTree_.matchedRecHits_electron_Dphi->at(idx).push_back(fabs(Dphi));
      outTree_.matchedRecHits_electron_Dz->at(idx).push_back(Dz);
      outTree_.matchedRecHits_electron_RDphi->at(idx).push_back(RDphi);
      outTree_.matchedRecHits_electron_DR->at(idx).push_back(DR);
      outTree_.matchedRecHits_electron_dist->at(idx).push_back(dist);
      
      if(recHit.energy() > 0.5)
      {
        sieie += energy*pow(eta-gp_track.eta(),2);
        sipip += energy*pow(phi-gp_track.phi(),2);
        ss_hit_count += energy;
      }
      
      ++recHitIt;
    }
    outTree_.matchedRecHits_sietaieta->at(idx).push_back( ss_hit_count>0 ? sqrt(sieie)/ss_hit_count : -999. );
    outTree_.matchedRecHits_siphiiphi->at(idx).push_back( ss_hit_count>0 ? sqrt(sipip)/ss_hit_count : -999. );
    
    
    //---find associated BTL clusters
    int clusterIt = 0;
    for(auto clusIt : clustersBTL)
    {    
      DetId id = clusIt.detId();
      const auto& det = mtdGeometry_ -> idToDet(id);
      const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
      const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());    
      for(auto cluster : clusIt)
      {
        MTDDetId mtdId(id);
        int RR = 0;
        int module = 0;
        int modType = 0;
        int crystal = 0;
        int ieta = 0;
        int iphi = 0;
        
        if( mtdId.mtdSubDetector() == MTDDetId::BTL )
        {
          BTLDetId btlId(id);
          RR = btlId.mtdRR();
          module = btlId.module();
          modType = btlId.modType();
          crystal = btlId.crystal();
          ieta = btlId.ieta(crysLayout_);
          iphi = btlId.iphi(crysLayout_);
        }
        
        double energy = cluster.energy();
        double time   = cluster.time();
        int size=cluster.size();
        int sizeX=cluster.sizeX();
        int sizeY=cluster.sizeY();
        
        MeasurementPoint mp(cluster.x(),cluster.y());
        LocalPoint lp = topo.localPosition(mp);
        GlobalPoint gp = det->toGlobal(lp);
        
        float eta = gp.eta();
        float phi = gp.phi();
        
        int closestPoint=-1;
        float minDist=999;
        for(unsigned int ic=0; ic<gp_ext.size();++ic)
        {
          GlobalPoint diff(gp_ext[ic].x()-gp.x(),gp_ext[ic].y()-gp.y(),gp_ext[ic].z()-gp.z());
          if (diff.mag()<minDist)
          {
            closestPoint=ic;
            minDist=diff.mag();
          }
        }
        
        if( closestPoint == -1 ) continue;
        
        GlobalPoint gp_track = gp_ext[closestPoint];      
        
        float Deta  = gp_track.mag() > 0. ? eta-gp_track.eta()           : -999.;
        float Dphi  = gp_track.mag() > 0. ? deltaPhi(phi,gp_track.phi()) : -999.;
        float DR    = gp_track.mag() > 0. ? sqrt(Deta*Deta+Dphi*Dphi)    : -999.;
        float Dz    = gp_track.mag() > 0. ? gp_track.z()-gp.z()          : -999.;
        float RDphi = gp_track.mag() > 0. ? sqrt(gp_track.perp2())*Dphi  : -999.;
        float dist  = gp_track.mag() > 0. ? (gp-gp_track).mag()          : -999.;
	
	
        if( DR > track_hit_DRMax_ || dist > track_hit_distMax_ ) continue;
        if( gp_track.mag() <= 0. ) continue;
        
        outTree_.matchedClusters_n->at(idx) += 1;
	
        outTree_.matchedClusters_idx->at(idx).push_back(idx);
        outTree_.matchedClusters_det->at(idx).push_back(1);
        outTree_.matchedClusters_energy->at(idx).push_back(energy);
        outTree_.matchedClusters_energyCorr->at(idx).push_back(energy*fabs(sin(trackRef->theta())));
        outTree_.matchedClusters_time->at(idx).push_back(time);
        outTree_.matchedClusters_rr->at(idx).push_back(RR);
        outTree_.matchedClusters_module->at(idx).push_back(module);
        outTree_.matchedClusters_modType->at(idx).push_back(modType);
        outTree_.matchedClusters_crystal->at(idx).push_back(crystal);
        outTree_.matchedClusters_ieta->at(idx).push_back(ieta);
        outTree_.matchedClusters_iphi->at(idx).push_back(iphi);
        outTree_.matchedClusters_size->at(idx).push_back(size);
        outTree_.matchedClusters_size_x->at(idx).push_back(sizeX);
        outTree_.matchedClusters_size_y->at(idx).push_back(sizeY);
        outTree_.matchedClusters_local_x->at(idx).push_back(lp.x());
        outTree_.matchedClusters_local_y->at(idx).push_back(lp.y());
        outTree_.matchedClusters_local_z->at(idx).push_back(lp.z());
        outTree_.matchedClusters_global_R->at(idx).push_back(sqrt(gp.perp2()));
        outTree_.matchedClusters_electron_Deta->at(idx).push_back(fabs(Deta));
        outTree_.matchedClusters_electron_Dphi->at(idx).push_back(fabs(Dphi));
        outTree_.matchedClusters_electron_Dz->at(idx).push_back(Dz);
        outTree_.matchedClusters_electron_RDphi->at(idx).push_back(RDphi);
        outTree_.matchedClusters_electron_DR->at(idx).push_back(DR);
        outTree_.matchedClusters_electron_dist->at(idx).push_back(dist);
	
        ++clusterIt;
      }
    }
    if( verbosity_ ) std::cout << "---\n\n\n" << std::endl;
    
    
    
    //---get associated ETL simHits
    if( verbosity_ ) std::cout << "*** simHits - n tot: " << simHitsETL.size() << std::endl;
    simHitIt = 0;
    for(auto simHit : simHitsETL)
    {
      ETLDetId id = simHit.detUnitId();
      DetId geoId = id.geographicalId();
      const auto& det = mtdGeometry_ -> idToDet(geoId);
      const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
      const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());
      
      double energy = simHit.energyLoss()*1000.;
      double time   = simHit.tof();
      
      if( (time < 0) || (time > 25) ) continue;
      
      int RR = id.mtdRR();
      int module = id.module();
      int modType = id.modType();
      int crystal = 0;
      int ieta = 0;
      int iphi = 0;
      
      // ETL is already in module-local coordinates so just scale to cm from mm
      Local3DPoint simscaled(0.1*simHit.entryPoint().x(),0.1*simHit.entryPoint().y(),0.1*simHit.entryPoint().z());
      const auto& thepixel = topo.pixel(simscaled); // mm -> cm here is the switch                               
      const uint8_t row(thepixel.first), col(thepixel.second);
      
      LocalPoint lp_entry(   simHit.entryPoint().x()/10.,   simHit.entryPoint().y()/10.,   simHit.entryPoint().z()/10.);
      LocalPoint lp_mid  (simHit.localPosition().x()/10.,simHit.localPosition().y()/10.,simHit.localPosition().z()/10.);
      LocalPoint lp_exit (    simHit.exitPoint().x()/10.,    simHit.exitPoint().y()/10.,    simHit.exitPoint().z()/10.);
      GlobalPoint gp_entry = det->toGlobal(topo.pixelToModuleLocalPoint(lp_entry ,row,col));
      GlobalPoint gp_mid   = det->toGlobal(topo.pixelToModuleLocalPoint(lp_mid   ,row,col));
      GlobalPoint gp_exit  = det->toGlobal(topo.pixelToModuleLocalPoint(lp_exit  ,row,col));
      
      float eta = gp_mid.eta();
      float phi = gp_mid.phi();

      int closestPoint=-1;
      float minDist=999;
      for(unsigned int ic=0; ic<gp_ext.size();++ic)
      {
        GlobalPoint diff(gp_ext[ic].x()-gp_mid.x(),gp_ext[ic].y()-gp_mid.y(),gp_ext[ic].z()-gp_mid.z());
        if (diff.mag()<minDist)
        {
          closestPoint=ic;
          minDist=diff.mag();
        }
      }
      
      if( closestPoint == -1 ) continue;
      
      GlobalPoint gp_track = gp_ext[closestPoint];      
      
      float Deta  = gp_track.mag() > 0. ? eta-gp_track.eta()           : -999.;
      float Dphi  = gp_track.mag() > 0. ? deltaPhi(phi,gp_track.phi()) : -999.;
      float DR    = gp_track.mag() > 0. ? sqrt(Deta*Deta+Dphi*Dphi)    : -999.;
      float Dz    = gp_track.mag() > 0. ? gp_track.z()-gp_mid.z()      : -999.;
      float RDphi = gp_track.mag() > 0. ? sqrt(gp_track.perp2())*Dphi  : -999.;
      float dist  = gp_track.mag() > 0. ? (gp_mid-gp_track).mag()      : -999.;
      if( DR < 0.05 && DR > 0. )
      {
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ":   energy: " << energy << " MeV   time: " << time << " ns"
                                   << "   RR: " << RR << "   module: " << module << "   modType: " << modType << "   crystal: " << crystal;
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ":   entryPoint:   local: " << PrintPosition(lp_entry)        << "   global: " << PrintPosition(gp_entry) << std::endl;
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ":     midPoint:   local: " << PrintPosition(lp_mid)          << "   global: " << PrintPosition(gp_mid)   << std::endl;
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ":    exitPoint:   local: " << PrintPosition(lp_exit)         << "   global: " << PrintPosition(gp_exit)  << std::endl;
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ":  detPosition:  global: " << PrintPosition(det->position()) << std::endl;
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ": hit position:   local: " << PrintPosition(lp_mid)          << "   global: " << PrintPosition(gp_mid) << "   DR: " << DR << "   dist: " << (gp_track-gp_mid).mag() << std::endl;
      }
      
      if( DR > track_hit_DRMax_ || dist > track_hit_distMax_ ) continue;
      if( gp_track.mag() <= 0. ) continue;
      
      outTree_.matchedSimHits_n->at(idx) += 1;
      
      outTree_.matchedSimHits_idx->at(idx).push_back(idx);
      outTree_.matchedSimHits_det->at(idx).push_back(2);
      outTree_.matchedSimHits_energy->at(idx).push_back(energy);
      outTree_.matchedSimHits_energyCorr->at(idx).push_back(energy*fabs(sin(trackRef->theta())));
      outTree_.matchedSimHits_time->at(idx).push_back(time);
      outTree_.matchedSimHits_rr->at(idx).push_back(RR);
      outTree_.matchedSimHits_module->at(idx).push_back(module);
      outTree_.matchedSimHits_modType->at(idx).push_back(modType);
      outTree_.matchedSimHits_crystal->at(idx).push_back(crystal);
      outTree_.matchedSimHits_ieta->at(idx).push_back(ieta);
      outTree_.matchedSimHits_iphi->at(idx).push_back(iphi);
      outTree_.matchedSimHits_entry_local_x->at(idx).push_back(lp_entry.x());
      outTree_.matchedSimHits_entry_local_y->at(idx).push_back(lp_entry.y());
      outTree_.matchedSimHits_entry_local_z->at(idx).push_back(lp_entry.z());
      outTree_.matchedSimHits_entry_global_R->at(idx).push_back(sqrt(gp_entry.perp2()));
      outTree_.matchedSimHits_exit_local_x->at(idx).push_back(lp_exit.x());
      outTree_.matchedSimHits_exit_local_y->at(idx).push_back(lp_exit.y());
      outTree_.matchedSimHits_exit_local_z->at(idx).push_back(lp_exit.z());
      outTree_.matchedSimHits_exit_global_R->at(idx).push_back(sqrt(gp_exit.perp2()));
      outTree_.matchedSimHits_electron_Deta->at(idx).push_back(fabs(Deta));
      outTree_.matchedSimHits_electron_Dphi->at(idx).push_back(fabs(Dphi));
      outTree_.matchedSimHits_electron_DR->at(idx).push_back(DR);
      outTree_.matchedSimHits_electron_Dz->at(idx).push_back(Dz);
      outTree_.matchedSimHits_electron_RDphi->at(idx).push_back(RDphi);
      outTree_.matchedSimHits_electron_dist->at(idx).push_back(dist);
      
      ++simHitIt;
    }
    if( verbosity_ ) std::cout << "---" << std::endl;
    
    
    //---find associated ETL recHits
    sieie=0; sipip=0;
    ss_hit_count=0;
    recHitIt = 0;
    if( verbosity_ ) std::cout << "*** recHits - tot: " << recHitsETL.size() << std::endl;
    for(auto recHit : recHitsETL)
    {
      ETLDetId id = recHit.id();
      DetId geoId = id.geographicalId(  );
      const auto& det = mtdGeometry_ -> idToDet(geoId);
      const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
      const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());    
      
      double energy = recHit.energy();
      double time   = recHit.time();
      
      MeasurementPoint mp(recHit.row()+0.5f,recHit.column()+0.5f);
      LocalPoint lp = topo.localPosition(mp);
      GlobalPoint gp = det->toGlobal(lp);

      float eta = gp.eta();
      float phi = gp.phi();
      
      int RR = id.mtdRR();
      int module = id.module();
      int modType = id.modType();
      int crystal = 0;
      int ieta = 0;
      int iphi = 0;
      
      int closestPoint=-1;
      float minDist=999;
      for (unsigned int ic=0; ic<gp_ext.size();++ic)
      {
        GlobalPoint diff(gp_ext[ic].x()-gp.x(),gp_ext[ic].y()-gp.y(),gp_ext[ic].z()-gp.z());
        if( diff.mag() < minDist )
        {
          closestPoint=ic;
          minDist=diff.mag();
        }
      }
      
      if( closestPoint == -1 ) continue;
      
      GlobalPoint gp_track = gp_ext[closestPoint];      
      
      float Deta  = gp_track.mag() > 0. ? eta-gp_track.eta()           : -999.;
      float Dphi  = gp_track.mag() > 0. ? deltaPhi(phi,gp_track.phi()) : -999.;
      float DR    = gp_track.mag() > 0. ? sqrt(Deta*Deta+Dphi*Dphi)    : -999.;
      float Dz    = gp_track.mag() > 0. ? gp_track.z()-gp.z()          : -999.;
      float RDphi = gp_track.mag() > 0. ? sqrt(gp_track.perp2())*Dphi  : -999.;
      float dist  = gp_track.mag() > 0. ? (gp-gp_track).mag()          : -999.;
      
      // if( DR < 0.2 && DR > 0. )
      {
        if( verbosity_ ) std::cout << ">>> " << recHitIt << ":   energy: " << energy << " MeV   time: " << time << " ns"
                                   << "   RR: " << RR << "   module: " << module << "   modType: " << modType << "   crystal: " << crystal;
        if( verbosity_ ) std::cout << ">>> " << recHitIt << ":  detPosition:  global: " << PrintPosition(det->position()) << "   DR: " << DR << "   dist: " << (gp_track-det->position()).mag() << std::endl;
        if( verbosity_ ) std::cout << ">>> " << recHitIt << ": hit position:   local: " << PrintPosition(lp) << "   global: " << PrintPosition(gp) << "DR: " << DR << "   dist: " << (gp_track-gp).mag() << std::endl;
        if( verbosity_ ) std::cout << ">>> " << recHitIt << ": extrapolated track " << abs(recHit.row()-int()) << " : " << PrintPosition(gp_track) << std::endl;
      }
      
      if( DR > track_hit_DRMax_ || dist > track_hit_distMax_ ) continue;
      if( gp_track.mag() <= 0. ) continue;
      
      outTree_.matchedRecHits_n->at(idx) += 1;
      
      outTree_.matchedRecHits_idx->at(idx).push_back(idx);
      outTree_.matchedRecHits_det->at(idx).push_back(2);
      outTree_.matchedRecHits_energy->at(idx).push_back(energy);
      outTree_.matchedRecHits_energyCorr->at(idx).push_back(energy*fabs(sin(trackRef->theta())));
      outTree_.matchedRecHits_time->at(idx).push_back(time);
      outTree_.matchedRecHits_rr->at(idx).push_back(RR);
      outTree_.matchedRecHits_module->at(idx).push_back(module);
      outTree_.matchedRecHits_modType->at(idx).push_back(modType);
      outTree_.matchedRecHits_crystal->at(idx).push_back(crystal);
      outTree_.matchedRecHits_ieta->at(idx).push_back(ieta);
      outTree_.matchedRecHits_iphi->at(idx).push_back(iphi);
      outTree_.matchedRecHits_local_x->at(idx).push_back(lp.x());
      outTree_.matchedRecHits_local_y->at(idx).push_back(lp.y());
      outTree_.matchedRecHits_local_z->at(idx).push_back(lp.z());
      outTree_.matchedRecHits_global_R->at(idx).push_back(sqrt(gp.perp2()));
      outTree_.matchedRecHits_electron_Deta->at(idx).push_back(fabs(Deta));
      outTree_.matchedRecHits_electron_Dphi->at(idx).push_back(fabs(Dphi));
      outTree_.matchedRecHits_electron_Dz->at(idx).push_back(Dz);
      outTree_.matchedRecHits_electron_RDphi->at(idx).push_back(RDphi);
      outTree_.matchedRecHits_electron_DR->at(idx).push_back(DR);
      outTree_.matchedRecHits_electron_dist->at(idx).push_back(dist);
      
      if(recHit.energy() > 0.5)
      {
        sieie += energy*pow(eta-gp_track.eta(),2);
        sipip += energy*pow(phi-gp_track.phi(),2);
        ss_hit_count += energy;
      }
      
      ++recHitIt;
    }
    outTree_.matchedRecHits_sietaieta->at(idx).push_back( ss_hit_count>0 ? sqrt(sieie)/ss_hit_count : -999. );
    outTree_.matchedRecHits_siphiiphi->at(idx).push_back( ss_hit_count>0 ? sqrt(sipip)/ss_hit_count : -999. );
    
    
    //---find associated ETL clusters
    clusterIt = 0;
    for(auto clusIt : clustersETL)
    {    
      DetId id = clusIt.detId();
      const auto& det = mtdGeometry_ -> idToDet(id);
      const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
      const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());    
      for ( auto cluster : clusIt)
      {
        MTDDetId mtdId(id);
        int RR = 0;
        int module = 0;
        int modType = 0;
        int crystal = 0;
        int ieta = 0;
        int iphi = 0;
	
        if( mtdId.mtdSubDetector() == MTDDetId::ETL )
        {
          ETLDetId btlId(id);
          RR = btlId.mtdRR();
          module = btlId.module();
          modType = btlId.modType();
          // crystal = btlId.crystal();
          // ieta = btlId.ieta();
          // iphi = btlId.iphi();
        }
	
        double energy = cluster.energy();
        double time   = cluster.time();
        int size=cluster.size();
        int sizeX=cluster.sizeX();
        int sizeY=cluster.sizeY();
	
        MeasurementPoint mp(cluster.x(),cluster.y());
        LocalPoint lp = topo.localPosition(mp);
        GlobalPoint gp = det->toGlobal(lp);
	
        float eta = gp.eta();
        float phi = gp.phi();
	
        int closestPoint=-1;
        float minDist=999;
        for(unsigned int ic=0; ic<gp_ext.size();++ic)
        {
          GlobalPoint diff(gp_ext[ic].x()-gp.x(),gp_ext[ic].y()-gp.y(),gp_ext[ic].z()-gp.z());
          if (diff.mag()<minDist)
          {
            closestPoint=ic;
            minDist=diff.mag();
          }
        }
	
        if( closestPoint == -1 ) continue;
        
        GlobalPoint gp_track = gp_ext[closestPoint];      
	
        float Deta  = gp_track.mag() > 0. ? eta-gp_track.eta()           : -999.;
        float Dphi  = gp_track.mag() > 0. ? deltaPhi(phi,gp_track.phi()) : -999.;
        float DR    = gp_track.mag() > 0. ? sqrt(Deta*Deta+Dphi*Dphi)    : -999.;
        float Dz    = gp_track.mag() > 0. ? gp_track.z()-gp.z()          : -999.;
        float RDphi = gp_track.mag() > 0. ? sqrt(gp_track.perp2())*Dphi  : -999.;
        float dist  = gp_track.mag() > 0. ? (gp-gp_track).mag()          : -999.;
	
        if( DR > track_hit_DRMax_ || dist > track_hit_distMax_ ) continue;
        if( gp_track.mag() <= 0. ) continue;
	
        outTree_.matchedClusters_n->at(idx) += 1;
	
        outTree_.matchedClusters_idx->at(idx).push_back(idx);
        outTree_.matchedClusters_det->at(idx).push_back(2);
        outTree_.matchedClusters_energy->at(idx).push_back(energy);
        outTree_.matchedClusters_energyCorr->at(idx).push_back(energy*fabs(sin(trackRef->theta())));
        outTree_.matchedClusters_time->at(idx).push_back(time);
        outTree_.matchedClusters_rr->at(idx).push_back(RR);
        outTree_.matchedClusters_module->at(idx).push_back(module);
        outTree_.matchedClusters_modType->at(idx).push_back(modType);
        outTree_.matchedClusters_crystal->at(idx).push_back(crystal);
        outTree_.matchedClusters_ieta->at(idx).push_back(ieta);
        outTree_.matchedClusters_iphi->at(idx).push_back(iphi);
        outTree_.matchedClusters_size->at(idx).push_back(size);
        outTree_.matchedClusters_size_x->at(idx).push_back(sizeX);
        outTree_.matchedClusters_size_y->at(idx).push_back(sizeY);
        outTree_.matchedClusters_local_x->at(idx).push_back(lp.x());
        outTree_.matchedClusters_local_y->at(idx).push_back(lp.y());
        outTree_.matchedClusters_local_z->at(idx).push_back(lp.z());
        outTree_.matchedClusters_global_R->at(idx).push_back(sqrt(gp.perp2()));
        outTree_.matchedClusters_electron_Deta->at(idx).push_back(fabs(Deta));
        outTree_.matchedClusters_electron_Dphi->at(idx).push_back(fabs(Dphi));
        outTree_.matchedClusters_electron_Dz->at(idx).push_back(Dz);
        outTree_.matchedClusters_electron_RDphi->at(idx).push_back(RDphi);
        outTree_.matchedClusters_electron_DR->at(idx).push_back(DR);
        outTree_.matchedClusters_electron_dist->at(idx).push_back(dist);
	
        ++clusterIt;
      }
    }
    if( verbosity_ ) std::cout << "---\n\n\n" << std::endl;
    
    
    ++idx;
  } //--- fill the tree - barrel electrons
  
  
  
  outTree_.GetTTreePtr()->Fill();
}


std::string FTLDumpElectrons::PrintPosition(const GlobalPoint& gp)
{
  std::stringstream output;
  
  output << "(";
  output << std::fixed << std::setprecision(3) << std::setw(8) << gp.x() << ",";
  output << std::fixed << std::setprecision(3) << std::setw(8) << gp.y() << ",";
  output << std::fixed << std::setprecision(3) << std::setw(8) << gp.z();
  output << ") cm";
  
  output << "   R: " << std::fixed << std::setprecision(3) << std::setw(7) << gp.perp();
  output << " cm";
  
  output << "   eta: " << std::setprecision(3) << std::setw(6) << gp.eta(); 
  output << "   phi: " << std::setprecision(3) << std::setw(6) << gp.phi();
  
  return output.str();
}

std::string FTLDumpElectrons::PrintPosition(const LocalPoint& lp)
{
  std::stringstream output;
  
  output << "(";
  output << std::fixed << std::setprecision(3) << std::setw(8) << lp.x() << ",";
  output << std::fixed << std::setprecision(3) << std::setw(8) << lp.y() << ",";
  output << std::fixed << std::setprecision(3) << std::setw(8) << lp.z();
  output << ") cm";
  
  return output.str();
}
DEFINE_FWK_MODULE(FTLDumpElectrons);

#endif
