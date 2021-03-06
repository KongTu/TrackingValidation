// -*- C++ -*-
//
// Package:    TrackingValidation/validation
// Class:      validation
// 
/**\class validation validation.cc TrackingValidation/validation/plugins/validation.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Zhoudunming Tu
//         Created:  Tue, 05 Apr 2016 08:12:26 GMT
//
//


// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <map>
#include <sstream>


#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TRandom.h>
#include <TNtuple.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

//////////////////////////////////////////////
// CMSSW user include files
#include "DataFormats/Common/interface/DetSetAlgorithm.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"


#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/CaloTowers/interface/CaloTower.h"  
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"  
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"   

// Particle Flow
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"

// Heavyion
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

// Particle Flow
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"

// Vertex significance
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"

// Root include files
#include "TTree.h"
//
// Track Matching and fake rate calculations     
//#include "RiceHIG/V0Analysis/interface/V0Validator.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class validation : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit validation(const edm::ParameterSet&);
      ~validation();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
  edm::EDGetTokenT<reco::Centrality> centralityToken_;
  edm::EDGetTokenT<int> centralityBinToken_;

  edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
  edm::EDGetTokenT<reco::TrackCollection> trackSrc_;
  edm::EDGetTokenT<reco::PFCandidateCollection> pfCandSrc_;
  edm::EDGetTokenT<CaloTowerCollection> towerSrc_;
  edm::EDGetTokenT<reco::GenParticleCollection> genSrc_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;

  edm::InputTag vertexName_;
  edm::InputTag trackName_;
  edm::InputTag pfCandName_;
  edm::InputTag towerName_;

  double offlineptErr_;
  double offlineDCA_;
  double offlineChi2_;
  double offlinenhits_;

  bool doCentrality_;
  bool doGenParticle_;

  TH1D* genEta;
  TH1D* genY;

  TH1D* nVtx;
  TH1D* vtxTracksSize;
  TH1D* vtxZ;
  TH1D* vtxX;
  TH1D* vtxY;

  TH1D* beamX;
  TH1D* beamY;
  TH1D* beamZ;

  TH1D* beamXerror;
  TH1D* beamYerror;
  TH1D* beamZerror;

  TH1D* pt;
  TH1D* ptError;
  TH1D* eta;
  TH1D* phi;
  TH1D* DCAz;
  TH1D* DCAxy;
  TH1D* Chi2n;
  TH1D* numberOfHits;
  TH1D* Algo;

  TH1D* Ntrk;
  TH1D* HFsum;
  TH1D* HFsumEt;
  TH1D* CaloEta;
  TH1D* cBins;
  TH2D* caloVsCbin;

  TH3D* pTvsEtaCent;
  TH3D* pTvsEtaSumET;


};



//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
validation::validation(const edm::ParameterSet& iConfig):
genSrc_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genSrc"))),
beamSpotSrc_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotSrc")))
{
   //now do what ever initialization is needed
  
  trackName_  =  iConfig.getParameter<edm::InputTag>("trackName");
  vertexName_ =  iConfig.getParameter<edm::InputTag>("vertexName");
  pfCandName_ =  iConfig.getUntrackedParameter<edm::InputTag>("pfCandName");
  towerName_ =  iConfig.getParameter<edm::InputTag>("towerName");

  trackSrc_ = consumes<reco::TrackCollection>(trackName_);
  vertexSrc_ = consumes<reco::VertexCollection>(vertexName_);
  pfCandSrc_ = consumes<reco::PFCandidateCollection>(pfCandName_);
  towerSrc_ = consumes<CaloTowerCollection>(towerName_);

  offlineptErr_ = iConfig.getUntrackedParameter<double>("offlineptErr", 0.0);
  offlineDCA_ = iConfig.getUntrackedParameter<double>("offlineDCA", 0.0);
  offlineChi2_ = iConfig.getUntrackedParameter<double>("offlineChi2", 0.0);
  offlinenhits_ = iConfig.getUntrackedParameter<double>("offlinenhits", 0.0);

  doCentrality_ = iConfig.getUntrackedParameter<bool>("doCentrality");
  doGenParticle_ = iConfig.getUntrackedParameter<bool>("doGenParticle");

  centralityToken_ = consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"));
  centralityBinToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("centralityBinSrc"));

}


validation::~validation()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
validation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;

  double hiBin_ = 0.;
  if( doCentrality_ ){
    edm::Handle<reco::Centrality> centrality;
    iEvent.getByToken(centralityToken_, centrality);
    edm::Handle<int> cbin;
    iEvent.getByToken(centralityBinToken_, cbin);
    
    hiBin_ = *cbin;
    cBins->Fill( hiBin_ );
  }

  if( doGenParticle_ ){
    edm::Handle<reco::GenParticleCollection> genParticleCollection;
    iEvent.getByToken(genSrc_, genParticleCollection);
    
    for(unsigned it=0; it<genParticleCollection->size(); ++it) {

      const reco::GenParticle & genCand = (*genParticleCollection)[it];
      int status = genCand.status();
      double geneta = genCand.eta();
      int gencharge = genCand.charge();
      double genrapidity = genCand.rapidity();

      if( status != 1 || gencharge == 0 ) continue;
      //if( fabs(geneta) > 2.4 ) continue;

      genEta->Fill( geneta );
      genY->Fill( genrapidity );

    }
  }

  edm::Handle<reco::BeamSpot> theBeamSpotHandle;
  iEvent.getByToken(beamSpotSrc_, theBeamSpotHandle);
  const reco::BeamSpot* theBeamSpot = theBeamSpotHandle.product();

  double bx = theBeamSpot->x0();
  double by = theBeamSpot->y0();
  double bz = theBeamSpot->z0();
  double bxError = theBeamSpot->x0Error();
  double byError = theBeamSpot->y0Error();
  double bzError = theBeamSpot->z0Error();

  beamX->Fill( bx );
  beamY->Fill( by );
  beamZ->Fill( bz );
  beamXerror->Fill( bxError );
  beamYerror->Fill( byError );
  beamZerror->Fill( bzError );


  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexSrc_,vertices);
  double bestvz=-999.9, bestvx=-999.9, bestvy=-999.9;
  double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
  const reco::Vertex & vtx = (*vertices)[0];
  bestvz = vtx.z(); 
  bestvx = vtx.x(); 
  bestvy = vtx.y();
  bestvzError = vtx.zError(); 
  bestvxError = vtx.xError(); 
  bestvyError = vtx.yError();
  
  //first selection; vertices
  if(bestvz < -15.0 || bestvz > 15.0) return;

  int num = 0;
  for( unsigned i = 0; i < vertices->size(); i++){
    const reco::Vertex & vtx1 = (*vertices)[i];
    if( vtx1.isFake() ) continue;
    vtxTracksSize->Fill( vtx1.tracksSize() );
    if( vtx1.tracksSize() < 2 ) continue;
    num++;

    
  }
  nVtx->Fill( num );

  if( vtx.tracksSize() < 2 || vtx.isFake() ) return;
  vtxZ->Fill( vtx.z() );
  vtxX->Fill( vtx.x() );
  vtxY->Fill( vtx.y() );

  Handle<CaloTowerCollection> towers;
  iEvent.getByToken(towerSrc_, towers);

  double energy = 0.;
  double sumEt = 0.;
  for(unsigned i = 0; i < towers->size(); ++i){

        const CaloTower & hit= (*towers)[i];
        if( fabs(hit.eta()) > 5.0 || fabs(hit.eta()) < 3.0  ) continue;
        CaloEta->Fill( hit.eta() );
        energy += (hit.hadEnergy() + hit.emEnergy());
        sumEt += (hit.hadEt()+hit.emEt());
  }
  
  HFsum->Fill( energy );
  HFsumEt->Fill( sumEt );

  Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(trackSrc_, tracks);
  if( !tracks.isValid() ) return;

  Handle<reco::PFCandidateCollection> pfCandidates;
  iEvent.getByToken(pfCandSrc_, pfCandidates);
  if( !pfCandidates.isValid() ) return;

  int total = 0;
  for(unsigned it = 0; it < tracks->size(); it++){

     const reco::Track & trk = (*tracks)[it];
  
     math::XYZPoint bestvtx(bestvx,bestvy,bestvz);

        double dzvtx = trk.dz(bestvtx);
        double dxyvtx = trk.dxy(bestvtx);
        double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
        double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError); 
        double nhits = trk.numberOfValidHits();
        double chi2n = trk.normalizedChi2();
        double nlayers = trk.hitPattern().trackerLayersWithMeasurement();
        chi2n = chi2n/nlayers;
        double algoOffline = trk.algo();

        if(!trk.quality(reco::TrackBase::highPurity)) continue;

        DCAz->Fill( fabs(dzvtx/dzerror) );
        DCAxy->Fill( fabs(dxyvtx/dxyerror) );
        ptError->Fill( fabs(trk.ptError())/trk.pt() );
        Chi2n->Fill( chi2n );
        numberOfHits->Fill( nhits );
        Algo->Fill( algoOffline );

        if(fabs(trk.ptError())/trk.pt() > offlineptErr_ ) continue;
        if(fabs(dzvtx/dzerror) > offlineDCA_) continue;
        if(fabs(dxyvtx/dxyerror) > offlineDCA_) continue;
        total++;//count multiplicity
        if(!(algoOffline==4 || algoOffline==6 || algoOffline==7 || algoOffline==5)) continue;
        if(chi2n > offlineChi2_) continue;
        if(nhits < offlinenhits_) continue;

        pTvsEtaCent->Fill(trk.eta(), trk.pt(), energy);
        pTvsEtaSumET->Fill(trk.eta(), trk.pt(), sumEt);

        pt->Fill( trk.pt() );
        eta->Fill( trk.eta() );
        phi->Fill( trk.phi() );

        double ecalEnergy = 0.;
        double hcalEnergy = 0.;
        for( unsigned ic = 0; ic < pfCandidates->size(); ic++ ) {//calo matching loops

          const reco::PFCandidate& cand = (*pfCandidates)[ic];

          int type = cand.particleId();

          // only charged hadrons and leptons can be asscociated with a track
          if(!(type == reco::PFCandidate::h ||     //type1
          type == reco::PFCandidate::e ||     //type2
          type == reco::PFCandidate::mu      //type3
          )) continue;

          reco::TrackRef trackRef = cand.trackRef();
          if( it == trackRef.key() ) {
            // cand_index = ic;
            ecalEnergy = cand.ecalEnergy();
            hcalEnergy = cand.hcalEnergy();
            break;
          }
        }

        double energy = ecalEnergy+hcalEnergy;
        double pT = energy/(TMath::CosH(trk.eta()));
        double ratio = pT/trk.pt();
        
        if( !doCentrality_ ) hiBin_ = 1.0;
        caloVsCbin->Fill(hiBin_, ratio);

  }


      
  Ntrk->Fill( total );

}// last braket 
// ------------ method called once each job just before starting event loop  ------------
void 
validation::beginJob()
{
  edm::Service<TFileService> fs;
    
  TH3D::SetDefaultSumw2();

  genEta = fs->make<TH1D>("genEta", ";#eta", 2000,-10,10);
  genY = fs->make<TH1D>("genY", ";y", 2000,-10,10);

  Ntrk = fs->make<TH1D>("Ntrk", ";Ntrk", 10000, 0,10000);
  HFsum = fs->make<TH1D>("HFsum", ";HFsum", 1000, 0, 10000);
  HFsumEt = fs->make<TH1D>("HFsumEt", ";HFsumEt", 1000, 0, 10000); 
  CaloEta = fs->make<TH1D>("CaloEta",";#eta", 1000,-5,5);
  cBins = fs->make<TH1D>("cBins",";cbins", 200, 0, 200);
  caloVsCbin = fs->make<TH2D>("caloVsCbin",";cbins;Et/pT",200,0,200,100,0,10);

  nVtx = fs->make<TH1D>("nVtx",";nVtx",15,0,15);
  vtxTracksSize = fs->make<TH1D>("vtxTracksSize",";vtxTracksSize",10000,0,10000);
  vtxZ = fs->make<TH1D>("vtxZ",";vtxZ",400,-20,20);
  vtxX = fs->make<TH1D>("vtxX",";vtxX",400,-2,2);
  vtxY = fs->make<TH1D>("vtxY",";vtxY",400,-2,2);

  beamX = fs->make<TH1D>("beamX", ";x0", 400,-20,20);
  beamY = fs->make<TH1D>("beamY", ";y0", 400,-20,20);
  beamZ = fs->make<TH1D>("beamZ", ";z0", 400,-20,20);

  beamXerror = fs->make<TH1D>("beamXerror", ";x0Error", 1000,0,1);
  beamYerror = fs->make<TH1D>("beamYerror", ";y0Error", 1000,0,1);
  beamZerror = fs->make<TH1D>("beamZerror", ";z0Error", 1000,0,1);

  pt = fs->make<TH1D>("pt",";pt",10000,0,1000);
  ptError = fs->make<TH1D>("ptError",";ptError",1000,0,1);
  eta = fs->make<TH1D>("eta",";eta",1000,-3.0,3.0);
  phi = fs->make<TH1D>("phi",";phi",1000,-4.0,4.0);
  DCAz = fs->make<TH1D>("DCAz",";DCAz",1000,0,100);
  DCAxy = fs->make<TH1D>("DCAxy",";DCAxy",1000,0,100);
  numberOfHits = fs->make<TH1D>("numberOfHits",";numberOfHits",30,0,30);
  Algo = fs->make<TH1D>("Algo",";Algo",20,0,20);
  Chi2n = fs->make<TH1D>("Chi2n",";Chi2n",1000,0,1);

  pTvsEtaCent = fs->make<TH3D>("pTvsEtaCent",";#eta;p_{T}(GeV);centrality", 6, -2.4,2.4, 1000,0,100,1000,0,10000);
  pTvsEtaSumET = fs->make<TH3D>("pTvsEtaSumET",";#eta;p_{T}(GeV);sumET", 6, -2.4,2.4, 1000,0,100,1000,0,10000);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
validation::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
validation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(validation);
