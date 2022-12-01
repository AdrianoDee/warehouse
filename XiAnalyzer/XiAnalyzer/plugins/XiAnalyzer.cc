// -*- C++ -*-
//
// Package:    XbFrame/XiAnalyzer
// Class:      XiAnalyzer
//
/**\class XiAnalyzer XiAnalyzer.cc XbFrame/XiAnalyzer/plugins/XiAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Sergey Polikarpov
//         Created:  Tue, 15 Aug 2017 01:04:12 GMT
//
//

// system include files
#include <memory>

/// framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/FWLite/interface/EventBase.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
/// triggers
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
// #include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
/// tracks
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

/// muons
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

/// vertex fits
#include "RecoVertex/VertexTools/interface/VertexDistance.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

//// gen ??
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/CLHEP/interface/Migration.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // for miniAOD
#include "DataFormats/TrackReco/interface/Track.h" // for miniAOD

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/Error.h"
#include "TFile.h"
#include "TTree.h"
#include "TVectorD.h"    // for fixing tracks
#include "TMatrixDSym.h" // for fixing tracks
#include "TH1.h"

#include <vector>
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TLorentzVector.h"
#include <utility>
#include <string>
#include <map>


using namespace edm;
using namespace std;
using namespace reco;
//
///////
///// UPDATED TO 2018 PDG PARTICLES !!
/// {{{
double PDG_MUON_MASS    =   0.1056583745;
double PDG_PION_MASS    =   0.13957061;
double PDG_PIOZ_MASS    =   0.1349770;
double PDG_KAON_MASS    =   0.493677;
double PDG_PROTON_MASS  =   0.9382720813;
double PDG_KSHORT_MASS  =   0.497611;
double PDG_KSHORT_DM    =   0.000013;
double PDG_KSHORT_TIME  =   0.8954 * 0.0000000001;
double PDG_KS_MASS      =   PDG_KSHORT_MASS;
double PDG_LAMBDA_MASS  =   1.115683 ;
double PDG_LAMBDA_DM    =   0.000006 ;
double PDG_LAMBDA_TIME  =   2.632 * 0.0000000001;
double PDG_SIGMA0_MASS  =   1.192642 ;
double PDG_XImunus_MASS =   1.32171 ;
double PDG_XImunus_DM   =   0.00007 ;
double PDG_XImunus_TIME =   1.639 * 0.0000000001;
double PDG_OMmunus_MASS =   1.67245 ;
double PDG_OMmunus_DM   =   0.00029 ;
double PDG_OMmunus_TIME =   0.821 * 0.0000000001;
double PDG_DPM_MASS     =   1.86965 ;
double PDG_DPM_DM       =   0.00005 ;
double PDG_DPM_TIME     =   1.040 * 0.000000000001 ;
double PDG_DZ_MASS      =   1.86483 ;
double PDG_DZ_DM        =   0.00005 ;
double PDG_DZ_TIME      =   0.4101 * 0.000000000001 ;
double PDG_DS_MASS      =   1.96834 ;
double PDG_DS_DM        =   0.00007 ;
double PDG_DS_TIME      =   0.504 * 0.000000000001 ;
double PDG_LAMCZ_MASS   =   2.28646 ;
double PDG_LAMCZ_DM     =   0.00031 ;
double PDG_LAMCZ_TIME   =   2.00 * 0.0000000000001;
double PDG_XICZ_MASS    =   2.47087 ;
double PDG_XICZ_DM      =   0.00031 ;
double PDG_XICZ_TIME    =   1.12 * 0.0000000000001;
double PDG_XICP_MASS    =   2.46787 ;
double PDG_XICP_DM      =   0.00030 ;
double PDG_XICP_TIME    =   4.42 * 0.0000000000001;
double PDG_KSTARZ_MASS  =   0.89555;
double PDG_KSTARZ_GAMMA =   0.0473 ;
double PDG_KSTARP_MASS  =   0.89176;
double PDG_KSTARP_GAMMA =   0.0503 ;
double PDG_PHI_MASS     =   1.019461;
double PDG_PHI_GAMMA    =   0.004249;
double PDG_JPSI_MASS    =   3.096900;
double PDG_PSI2S_MASS   =   3.686097;
double PDG_X3872_MASS   =   3.87169;
double PDG_BU_MASS      =   5.27932;
double PDG_BU_TIME      =   1.638 * 0.000000000001;
double PDG_B0_MASS      =   5.27963;
double PDG_B0_TIME      =   1.520 * 0.000000000001;
double PDG_BS_MASS      =   5.36689;
double PDG_BS_TIME      =   1.509 * 0.000000000001;
double PDG_BC_MASS      =   6.2749;
double PDG_BC_TIME      =   0.507 * 0.000000000001;
double PDG_LB_MASS      =   5.61960;
double PDG_LB_TIME      =   1.470 * 0.000000000001;
double PDG_XIBZ_MASS    =   5.7919;
double PDG_XIBZ_TIME    =   1.479 * 0.000000000001;
double PDG_XIBM_MASS    =   5.7970;
double PDG_XIBM_TIME    =   1.571 * 0.000000000001;
double PDG_OMBM_MASS    =   6.0461;
double PDG_OMBM_TIME    =   1.64 * 0.000000000001;
double PDG_C            =   29979245800.; // in cm/c
//
/// }}}

ParticleMass PM_PDG_MUON_MASS = PDG_MUON_MASS;
ParticleMass PM_PDG_JPSI_MASS = PDG_JPSI_MASS;
ParticleMass PM_PDG_KAON_MASS = PDG_KAON_MASS;
ParticleMass PM_PDG_PION_MASS = PDG_PION_MASS;
ParticleMass PM_PDG_PROTON_MASS = PDG_PROTON_MASS;
ParticleMass PM_PDG_PSI2S_MASS = PDG_PSI2S_MASS;

reco::Track fix_track(const reco::Track *tk, double delta=1e-8);

reco::Track fix_track(const reco::TrackRef& tk)
{
    reco::Track t = reco::Track(*tk);
    return fix_track(&t);
}

/* Check for a not positive definite covariance matrix. If the covariance matrix is not positive definite, we force it to be positive definite by
 * adding the minimum eigenvalue to the diagonal of the covariance matrix plus `delta`.
 * See https://nhigham.com/2020/12/22/what-is-a-modified-cholesky-factorization/ */
reco::Track fix_track(const reco::Track *tk, double delta)
{
    unsigned int i, j;
    double min_eig = 1;

    /* Get the original covariance matrix. */
    reco::TrackBase::CovarianceMatrix cov = tk->covariance();

    /* Convert it from an SMatrix to a TMatrixD so we can get the eigenvalues. */
    TMatrixDSym new_cov(cov.kRows);
    for (i = 0; i < cov.kRows; i++) {
        for (j = 0; j < cov.kRows; j++) {
            /* Need to check for nan or inf, because for some reason these
             * cause a segfault when calling Eigenvectors().
             *
             * No idea what to do here or why this happens. */
            if (std::isnan(cov(i,j)) || std::isinf(cov(i,j)))
                cov(i,j) = 1e-6;
            new_cov(i,j) = cov(i,j);
        }
    }

    /* Get the eigenvalues. */
    TVectorD eig(cov.kRows);
    new_cov.EigenVectors(eig);
    for (i = 0; i < cov.kRows; i++)
        if (eig(i) < min_eig)
            min_eig = eig(i);

    /* If the minimum eigenvalue is less than zero, then subtract it from the
     * diagonal and add `delta`. */
    if (min_eig < 0) {
        for (i = 0; i < cov.kRows; i++)
            cov(i,i) -= min_eig - delta;
    }

    return reco::Track(tk->chi2(), tk->ndof(), tk->referencePoint(), tk->momentum(), tk->charge(), cov, tk->algo(), (reco::TrackBase::TrackQuality) tk->qualityMask());
}

//
//
// class declaration
//


class XiAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit XiAnalyzer(const edm::ParameterSet&);
      ~XiAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
// HLTConfigProvider hltConfig_;
edm::EDGetTokenT<reco::BeamSpot> thebeamspot_;
edm::EDGetTokenT<reco::VertexCollection> vtxSample;

// edm::EDGetTokenT<vector < pat::GenericParticle > > tracks_; /// AOD OLNY
edm::EDGetTokenT<pat::PackedCandidateCollection> trkTkn; /// shiAOD
edm::EDGetTokenT<pat::PackedCandidateCollection> trkTkndisc; /// shiAOD
edm::EDGetTokenT<pat::PackedCandidateCollection> trkTknlost; /// shiAOD
edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> tok_v0_; /// shiAOD
edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;

edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> builderToken_;

std::vector<float> *LA_mass     , *LA_vtxp;
std::vector<int>   *LA_pr_charg;
std::vector<float> *LA_px_CL    , *LA_py_CL     , *LA_pz_CL;
std::vector<float> *LA_VtxX_CL  , *LA_VtxY_CL   , *LA_VtxZ_CL;
std::vector<float> *LA_VtxXE_CL , *LA_VtxYE_CL  , *LA_VtxZE_CL;
std::vector<float> *LA_pi_px    , *LA_pi_py     , *LA_pi_pz;
std::vector<float> *LA_pr_px    , *LA_pr_py     , *LA_pr_pz;
std::vector<float> *LA_pi_ips;

std::vector<float> *XI_mass     , *XI_vtxp;
std::vector<int>   *XI_pi_charg , *XI_pi_purit  , *XI_pi_lost;
std::vector<float> *XI_px       , *XI_py        , *XI_pz;
std::vector<float> *XI_pi_px    , *XI_pi_py     , *XI_pi_pz;
std::vector<float> *XI_pi_px_CV , *XI_pi_py_CV  , *XI_pi_pz_CV;
std::vector<float> *XI_VtxX     , *XI_VtxY      , *XI_VtxZ;
std::vector<float> *XI_VtxXE    , *XI_VtxYE     , *XI_VtxZE;
std::vector<int>   *XI_pi_hit   , *XI_pi_pix;
std::vector<float> *XI_pi_ips;

std::vector<float> *PV_becos_XX , *PV_becos_YY  , *PV_becos_ZZ;
std::vector<float> *PV_becos_EX , *PV_becos_EY  , *PV_becos_EZ;
std::vector<float> *PV_becos_CL;
std::vector<int>   *PV_becos_dN;

// std::vector<TLorentzVector> *gen_xim_p4;
std::vector<float> *gen_XI_px   , *gen_XI_py    , *gen_XI_pz, *gen_XI_ma;

Int_t       nCand;

Int_t       run;
Int_t       event;
Float_t     lumi;

Int_t       numPV;
Int_t       numTrack;
Int_t       numLost;
Int_t       numLosa;
Int_t       numLosx;
Int_t       numV0;

TTree *wwtree;
// TTree *wwtree1;
TFile *f;
std::string fileName;

TLorentzVector gen_xim_p4;

};

//
XiAnalyzer::XiAnalyzer(const edm::ParameterSet& iConfig) :

LA_mass(0)      , LA_vtxp(0)    , LA_pr_charg(0),
LA_px_CL(0)     , LA_py_CL(0)   , LA_pz_CL(0),
LA_VtxX_CL(0)   , LA_VtxY_CL(0) , LA_VtxZ_CL(0),
LA_VtxXE_CL(0)  , LA_VtxYE_CL(0), LA_VtxZE_CL(0),
LA_pi_px(0)     , LA_pi_py(0)   , LA_pi_pz(0),
LA_pr_px(0)     , LA_pr_py(0)   , LA_pr_pz(0),
LA_pi_ips(0)    ,

XI_mass(0)      , XI_vtxp(0)    , XI_pi_charg(0), XI_pi_purit(0), XI_pi_lost(0),
XI_px(0)        , XI_py(0)      , XI_pz(0),
XI_pi_px(0)     , XI_pi_py(0)   , XI_pi_pz(0),
XI_pi_px_CV(0)  , XI_pi_py_CV(0), XI_pi_pz_CV(0),
XI_VtxX(0)      , XI_VtxY(0)    , XI_VtxZ(0),
XI_VtxXE(0)     , XI_VtxYE(0)   , XI_VtxZE(0),
XI_pi_hit(0)    , XI_pi_pix(0)  , XI_pi_ips(0),

PV_becos_XX(0)  , PV_becos_YY(0), PV_becos_ZZ(0),
PV_becos_EX(0)  , PV_becos_EY(0), PV_becos_EZ(0),
PV_becos_CL(0)  , PV_becos_dN(0),

gen_XI_px(0)    , gen_XI_py(0)  , gen_XI_pz(0), gen_XI_ma(0),
//
nCand(0),

run(0),
event(0),
lumi(0),

numPV(0),
numTrack(0),
numLost(0),
numLosa(0),
numLosx(0),
numV0(0)
{
// fileName = iConfig.getUntrackedParameter<std::string>("fileName","BFinder.root");
fileName = iConfig.getUntrackedParameter<std::string>("fileName");
// fileName = iConfig.getUntrackedParameter<std::string>("fileName","BFinder_vR.root");
   //now do what ever initialization is needed
   usesResource("TFileService");
// tok_v0_     = consumes<reco::VertexCompositeCandidateCollection>(   edm::InputTag("generalV0Candidates:Lambda")); /// AOD
vtxSample       =   consumes<reco::VertexCollection>(edm::InputTag("offlineSlimmedPrimaryVertices")); /// shiAOD
thebeamspot_    =   consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));  /// shiAOD
trkTkn          =   consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidates")); // shiAOD
trkTkndisc      =   consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidatesDiscarded")); /// shiAOD
trkTknlost      =   consumes<pat::PackedCandidateCollection>(edm::InputTag("lostTracks")); ///
tok_v0_         =   consumes<reco::VertexCompositePtrCandidateCollection>(edm::InputTag("slimmedLambdaVertices")); /// shiAOD
// packedGenToken_ =   consumes<edm::View<pat::PackedGenParticle> >(edm::InputTag("packedGenParticles"));
builderToken_ = esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"));
}


XiAnalyzer::~XiAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
XiAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
using namespace edm;
using namespace reco;
using namespace std;
using reco::MuonCollection;

run   = iEvent.id().run();
event = iEvent.id().event();

lumi = 1;
lumi = iEvent.luminosityBlock();


// edm::Handle<reco::GenParticleCollection> GenParticles;
// iEvent.getByToken(genToken_, GenParticles);
// gen_xim_p4.SetPtEtaPhiM(0.,0.,0.,0.);

/*
Handle<edm::View<pat::PackedGenParticle> > packed;
iEvent.getByToken(packedGenToken_,packed);
//

// std::cout << run << ":" << event << " ALARM !!! "<< packed->size() << std::endl;
int _numgen = packed->size() ;
int eventpassed = 0;
//
if ( packed.isValid() )
{
    for(Int_t ii=0; ii<_numgen; ii++)
    {
        Int_t pdgId = (*packed)[ii].pdgId();
	    if (abs(pdgId) == 3312)
        {
            // std::cout << run << ":" << event << " ALARM !!! "<< (*packed)[ii].pt() << std::endl;
            if ((*packed)[ii].pt() < 1.0) continue;
            if (fabs((*packed)[ii].eta()) > 3.0) continue;
            //
            TLorentzVector p4xi_gen;
	        p4xi_gen.SetPtEtaPhiM((*packed)[ii].pt(),(*packed)[ii].eta(),(*packed)[ii].phi(),(*packed)[ii].mass());
	        /// gen_xim_p4->push_back(p4xi_gen);
	        gen_XI_px->push_back    (p4xi_gen.X());
	        gen_XI_py->push_back    (p4xi_gen.Y());
	        gen_XI_pz->push_back    (p4xi_gen.Z());
	        gen_XI_ma->push_back    (p4xi_gen.M());
            //
            eventpassed += 1;
        }
    }
}

if (eventpassed == 0)
{
    // std::cout << run << ":" << event << " ALARM !!! EVENT NOT PASSED " << std::endl << std::endl;
    //return;
}

*/
// if (eventpassed > 1)
// {
//     std::cout << run << ":" << event << " ALARM !!! more than one GEN XI " << std::endl << std::endl;
//     //return;
// }
// std::cout << " fff 1111 fff" << std::endl;

// wwtree1->Fill();

// declare new track builder for new Transient track collection for miniAOD ???


// ESHandle<TransientTrackBuilder> theB;
// iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

auto const &theB = iSetup.getData(builderToken_);

// Get HLT results HLT   HLT  HLT   HLT  HLT   HLT
// edm::Handle<edm::TriggerResults> triggerResults_handle;
// iEvent.getByToken(hlTriggerResults_, triggerResults_handle);

//// PV PV  PV  PV  PV  PV  PV  PV  PV  PV  PV  PV  PV
Handle < VertexCollection > recVtxs;
iEvent.getByToken(vtxSample, recVtxs);
///
Vertex thePrimaryV;
thePrimaryV = Vertex(*(recVtxs->begin()));
const reco::VertexCollection & vertices = *recVtxs.product();
// for (reco::VertexCollection::const_iterator vtx = recVtxs->begin(); vtx != recVtxs->end(); ++vtx)
// {
//     nVtxTrks = vtx->tracksSize();
//     thePrimaryV = Vertex(*vtx);
// }


///  MUONS + TRACKS    MUONS + TRACKS    MUONS + TRACKS
// Handle < vector < pat::GenericParticle > >thePATTrackHandle; /// AOD ONLY
// iEvent.getByToken(tracks_, thePATTrackHandle); /// AOD ONLY
Handle<pat::PackedCandidateCollection> tracks; /// shiAOD
Handle<pat::PackedCandidateCollection> discTracks; /// shiAOD
Handle<pat::PackedCandidateCollection> lostTracks; /// shiAOD
iEvent.getByToken(trkTkn, tracks); /// shiAOD
iEvent.getByToken(trkTkndisc, discTracks); /// shiAOD
iEvent.getByToken(trkTknlost, lostTracks); /// shiAOD

numTrack    = tracks->size();
numLost     = lostTracks->size();

////  V0  V0  V0  V0  V0  V0  V0  V0  V0  V0  V0  V0
Handle<reco::VertexCompositePtrCandidateCollection> v0Coll; // shiAOD
iEvent.getByToken(tok_v0_, v0Coll);

numPV       = vertices.size();
numV0       = v0Coll->size();

// if (numV0   < 1 ) return;
// if (numTrack< 1 ) return;

numLosa     = 0;
numLosx     = 0;
std::set<int> trkids;
std::set<int> lotids;
std::set<int> loXids;

// cout << "numV0 " << numV0 << " numTrack " << numTrack << "\n";

// int trig = 0; // bit-based variable for trigger description

///// 2017 , 2018. not 2016 !
// unsigned int NTRIGGERS = 14;
// std::string TriggersToTest[NTRIGGERS] = {
// "HLT_Dimuon25_Jpsi","HLT_Dimuon20_Jpsi_Barrel_Seagulls", // 0, 1: inclusive dimuon jpsi
// "HLT_DoubleMu4_JpsiTrk_Displaced","HLT_DoubleMu4_JpsiTrkTrk_Displaced", // 2, 3: displaced jpistrk or jpsi trktrk
// "HLT_DoubleMu4_3_Jpsi_Displaced", "HLT_DoubleMu4_3_Jpsi", "HLT_DoubleMu4_Jpsi_Displaced",  // 4, 5, 6: prescaled, 4 for 2017, 5 for 2018
// "HLT_Dimuon18_PsiPrime", "HLT_Dimuon10_PsiPrime_Barrel_Seagulls", // 7, 8: inclusive dimuon psi2s
// "HLT_DoubleMu4_PsiPrimeTrk_Displaced", // 9: displaced psi2s trk
// "HLT_Dimuon0_Jpsi3p5_Muon2", // 10: triple-mu (jpsi + muon)
// "HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi", // 11: jpsi + 2 trkmu (phi->mumu)
// "HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi", "HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05" // 12, 13: jpsi+2 trk(phi->KK), 12 for 2017, 13 for 2018
// };
//
// std::string TriggerFilters[NTRIGGERS] = {
// "hltDisplacedmumuFilterDimuon25Jpsis","hltDisplacedmumuFilterDimuon20JpsiBarrelnoCow", // 0, 1: inclusive dimuon jpsi
// "hltJpsiTkVertexFilter","hltJpsiTkTkVertexFilterPhiKstar", // 2, 3: displaced jpistrk or jpsi trktrk
// "hltDisplacedmumuFilterDoubleMu43Jpsi", "hltmumuFilterDoubleMu43Jpsi", "hltDisplacedmumuFilterDoubleMu4Jpsi",  // 4, 5, 6: prescaled, 4 for 2017, 5 for 2018
// "hltDisplacedmumuFilterDimuon18PsiPrimes", "hltDisplacedmumuFilterDimuon10PsiPrimeBarrelnoCow", // 7, 8: inclusive dimuon psi2s
// "hltPsiPrimeTkVertexFilter", // 9: displaced psi2s trk
// "hltVertexmumuFilterJpsiMuon3p5", // 10: triple-mu (jpsi + muon)
// "hltDiMuonGlbOrTrk0zFiltered0p2v2", // 11: jpsi + 2 trkmu (phi->mumu)
// "hltJpsiTkTkVertexFilterPhiDoubleTrk1v2", "hltJpsiTkTkVertexFilterPhiDoubleTrk1v4" // 12, 13: jpsi+2 trk(phi->KK), 12 for 2017, 13 for 2018
// };

// unsigned short int TriggersFired[NTRIGGERS] = {};
//
// std::string strfired;
//
// if (triggerResults_handle.isValid())
// {
//     const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
//     for (unsigned int i = 0; i < NTRIGGERS; i++)
//     {
//         for (int version = 1; version < 30; version++)
//         {
//             std::stringstream ss; // full trigger name
//             ss << TriggersToTest[i] << "_v" << version;
//             //
//             unsigned int bit = TheTriggerNames.triggerIndex(edm::InputTag(ss.str()).label());
//             if (bit < triggerResults_handle->size() && triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit))
//             {
//                 trig += (1<<i);
//                 TriggersFired[i] = 1;
//                 strfired.append( ss.str());
//                 strfired.append( " " );
//                 //break;
//             }
//         }
//     }
// }
// else
// {
//     std::cout << " No trigger Results in event :( " << run << "," << event << std::endl;
// }
//

nCand = 0;

float PM_sigma = 1.e-7;
float chi = 0.;
float ndf = 0.;
int fitgood = 1;
KinematicParticleFactoryFromTransientTrack pFactory;

for (reco::VertexCompositePtrCandidateCollection::const_iterator vee = v0Coll->begin(); vee != v0Coll->end(); vee++) {
    TransientTrack pionTT;
    TransientTrack protTT;
    //check momentum and add pion first, prot second
    // std::cout << " starting a lambda 0 " << run << "," << event <<  " "  << std::endl;
    //
    fitgood = 1;
    try
    {
        if ((vee->daughter(0)->bestTrack() == nullptr) || (vee->daughter(1)->bestTrack() == nullptr))
        {
            fitgood=0;
        }
        else
        {
            if (vee->daughter(0)->momentum().mag2() < vee->daughter(1)->momentum().mag2())
            {
                pionTT = theB.build( fix_track(vee->daughter(0)->bestTrack()) );
                protTT = theB.build( fix_track(vee->daughter(1)->bestTrack()) );
            }
            else
            {
                pionTT = theB.build( fix_track(vee->daughter(1)->bestTrack()) );
                protTT = theB.build( fix_track(vee->daughter(0)->bestTrack()) );
            }
        }
    }
    catch (const VertexException &)
    {
        fitgood = 0;
    }
    // std::cout << " starting a lambda 1 " << run << "," << event <<  " "  << std::endl;
    if (fitgood == 0) continue;
    // std::cout << " starting a lambda 1.5 " << run << "," << event <<  " "  << std::endl;
    if(!pionTT.isValid()) continue;
    if(!protTT.isValid()) continue;
    // std::cout << " starting a lambda 2 " << run << "," << event <<  " "  << std::endl;

    vector<RefCountedKinematicParticle> VeeParticles;
    VeeParticles.push_back(pFactory.particle(pionTT,  PM_PDG_PION_MASS,     chi,ndf,PM_sigma));
    VeeParticles.push_back(pFactory.particle(protTT,  PM_PDG_PROTON_MASS,   chi,ndf,PM_sigma));

    // std::cout << " starting a lambda 3 " << run << "," << event <<  " "  << std::endl;
    KinematicParticleVertexFitter Lfitter;
    RefCountedKinematicTree veeVertexFitTree;
    fitgood = 1;
    try
    {
        veeVertexFitTree = Lfitter.fit(VeeParticles); // v0 without constraints
    }
    catch (const VertexException &)
    {
        fitgood = 0;
    }
    if (fitgood == 0) continue;
    if (!veeVertexFitTree->isValid())  {continue; }

    veeVertexFitTree->movePointerToTheTop();
    RefCountedKinematicParticle vee_vFit_noC  = veeVertexFitTree->currentParticle();
    RefCountedKinematicVertex veeVtx_vFit_noC = veeVertexFitTree->currentDecayVertex();
    double MY_lamass2 = vee_vFit_noC->currentState().mass(); // lamb without constraints;

    veeVertexFitTree->movePointerToTheFirstChild();
    RefCountedKinematicParticle pionCand    = veeVertexFitTree->currentParticle();
    veeVertexFitTree->movePointerToTheNextChild();  // her daughters
    RefCountedKinematicParticle protCand    = veeVertexFitTree->currentParticle();
    veeVertexFitTree->movePointerToTheTop();
//
    ParticleMass lamb_mass = PDG_LAMBDA_MASS;  float lamb_sigma = PDG_LAMBDA_DM;
    KinematicParticleFitter csFitterVee;
    KinematicConstraint * lam_c = new MassKinematicConstraint(lamb_mass, lamb_sigma);
    // add mass constraint to the lamb fit to do a constrained fit:

    veeVertexFitTree = csFitterVee.fit(lam_c,veeVertexFitTree);
    if (!veeVertexFitTree->isValid())  { continue;  }

    veeVertexFitTree->movePointerToTheTop();
    RefCountedKinematicParticle vee_vFit_withMC = veeVertexFitTree->currentParticle();
    RefCountedKinematicVertex   veeVtx_vFit_withMC = veeVertexFitTree->currentDecayVertex();
    double LA_Prob_tmp   = TMath::Prob(veeVtx_vFit_withMC->chiSquared(), (int) veeVtx_vFit_withMC->degreesOfFreedom());
    if (LA_Prob_tmp < 0.01) continue;
    // std::cout << "mass 2 " <<  MY_lamass2 << std::endl;
    //
    TLorentzVector p4la;
    p4la.SetXYZM(vee_vFit_withMC->currentState().globalMomentum().x(), vee_vFit_withMC->currentState().globalMomentum().y(), vee_vFit_withMC->currentState().globalMomentum().z(), vee_vFit_withMC->currentState().mass() );
    //
    //

      LA_mass         ->push_back( MY_lamass2                                             );
        LA_vtxp         ->push_back( LA_Prob_tmp                                            );
        LA_pr_charg     ->push_back( protTT.charge()                                        );
        LA_px_CL        ->push_back( vee_vFit_withMC->currentState().globalMomentum().x()   );
        LA_py_CL        ->push_back( vee_vFit_withMC->currentState().globalMomentum().y()   );
        LA_pz_CL        ->push_back( vee_vFit_withMC->currentState().globalMomentum().z()   );
        LA_VtxX_CL      ->push_back( veeVtx_vFit_withMC->position().x()                     );
        LA_VtxY_CL      ->push_back( veeVtx_vFit_withMC->position().y()                     );
        LA_VtxZ_CL      ->push_back( veeVtx_vFit_withMC->position().z()                     );
        LA_VtxXE_CL     ->push_back( veeVtx_vFit_withMC->error().cxx()                      );
        LA_VtxYE_CL     ->push_back( veeVtx_vFit_withMC->error().cyy()                      );
        LA_VtxZE_CL     ->push_back( veeVtx_vFit_withMC->error().czz()                      );
        LA_pi_px        ->push_back( pionCand->currentState().globalMomentum().x()          );
        LA_pi_py        ->push_back( pionCand->currentState().globalMomentum().y()          );
        LA_pi_pz        ->push_back( pionCand->currentState().globalMomentum().z()          );
        LA_pr_px        ->push_back( protCand->currentState().globalMomentum().x()          );
        LA_pr_py        ->push_back( protCand->currentState().globalMomentum().y()          );
        LA_pr_pz        ->push_back( protCand->currentState().globalMomentum().z()          );
        //LA_pi_ips       ->push_back( fabs(vee->daughter(0)->bestTrack()->dxy(bestVtxBSIP.position())) / (0.000001 + fabs(vee->daughter(0)->bestTrack()->dxyError())) );

continue;
    double Xi_Prob_tmp   = 0.0000001;
    //
    int trknumb = 0;
    /// LOOPING OVER NORMAL TRACKS WOHOOO
    for (pat::PackedCandidateCollection::const_iterator iTra1 = tracks->begin(); iTra1 != tracks->end(); iTra1++)
    {
        if (!(iTra1->hasTrackDetails())) continue;
        auto iTrack1 = iTra1->bestTrack();
        if (iTrack1 == nullptr) continue;
        //
        trkids.insert(trknumb);
        trknumb += 1;
        //
        if(iTrack1->pt()<0.2) continue;
        if( fabs(iTrack1->eta()) > 3.0 ) continue;
        //
        if(iTrack1 == vee->daughter(0)->bestTrack() ) continue;
        if(iTrack1 == vee->daughter(0)->bestTrack() ) continue;
        if (iTrack1->charge() * protTT.charge() > 0) continue;
        // if (!(iTrack1->quality(reco::TrackBase::highPurity))) continue;
        int pi1quality = 0; // 1 = loose, 2 = tight, 4 = highPurity
        if(iTrack1->quality(reco::TrackBase::loose))
        {
            pi1quality += 1;
        }
        if(iTrack1->quality(reco::TrackBase::tight))
        {
            pi1quality += 2;
        }
        if(iTrack1->quality(reco::TrackBase::highPurity))
        {
            pi1quality += 4;
        }
        //
        TransientTrack pion1TT = theB.build( fix_track(iTrack1) );
        if (!pion1TT.isValid()) continue;
        if (pionTT == pion1TT) continue;
        TLorentzVector p4pi1;
        p4pi1.SetPtEtaPhiM(iTrack1->pt(),iTrack1->eta(),iTrack1->phi(), PDG_PION_MASS);
        //
        if ((p4pi1 + p4la).M() > 2.0) continue; // D0 selection, + pimass = +0.135 = 2.0
        //
        //
        std::vector<RefCountedKinematicParticle> Xi_candidate_init;
        Xi_candidate_init.push_back(pFactory.particle(pion1TT, PM_PDG_PION_MASS, chi,ndf, PM_sigma));
        Xi_candidate_init.push_back(vee_vFit_withMC);
        RefCountedKinematicTree XiTree;
        KinematicParticleVertexFitter xFitter; //KinematicParticleVertexFitter
        fitgood = 1;
        try
        {
            XiTree = xFitter.fit(Xi_candidate_init);
        }
        catch (const VertexException &)
        {
            fitgood = 0;
        }
        if (fitgood == 0) continue;
        //
        if (!XiTree->isValid()) continue;
        XiTree->movePointerToTheTop();
        RefCountedKinematicParticle Xi_vFit = XiTree->currentParticle();
        RefCountedKinematicVertex   Xi_vtx = XiTree->currentDecayVertex();
        //
        double Xi_mass_tmp = Xi_vFit->currentState().mass();
        if (Xi_mass_tmp > 1.36) continue; //
        if (Xi_vtx->chiSquared() < 0) continue;
        Xi_Prob_tmp   = (double) (0.0 + TMath::Prob(Xi_vtx->chiSquared(), Xi_vtx->degreesOfFreedom()));
        if(Xi_Prob_tmp < 0.01) continue;
        //
        XiTree->movePointerToTheFirstChild();
        RefCountedKinematicParticle pixiCand    = XiTree->currentParticle();
        XiTree->movePointerToTheTop();
        //
        //
        // topological requirements
        double tmp_dx   = veeVtx_vFit_withMC->position().x() - Xi_vtx->position().x() ;
        double tmp_dy   = veeVtx_vFit_withMC->position().y() - Xi_vtx->position().y() ;
        double tmp_edx  = fabs(veeVtx_vFit_withMC->error().cxx()) + fabs(Xi_vtx->error().cxx()) ;
        double tmp_edy  = fabs(veeVtx_vFit_withMC->error().cyy()) + fabs(Xi_vtx->error().cyy()) ;
        double tmp_px   = vee_vFit_withMC->currentState().globalMomentum().x();
        double tmp_py   = vee_vFit_withMC->currentState().globalMomentum().y();
        //
        double tmp_DS2  = sqrt((tmp_dx*tmp_dx)/tmp_edx + (tmp_dy*tmp_dy)/tmp_edy );
        if( tmp_DS2 < 2) continue;
        double tmp_Cos2 = (tmp_dx * tmp_px + tmp_dy * tmp_py) / (sqrt(tmp_dx * tmp_dx + tmp_dy * tmp_dy)*sqrt(tmp_px * tmp_px + tmp_py * tmp_py ));
        if( tmp_Cos2 < 0.5) continue;
        //
        ////
        reco::Vertex bestVtxBSIP;  /// using Xi_vFit and Xi_vtx HERE !!
         // {{{ GET THE BEST PV BY CHOSING THE BEST POINTING ANGLE AND REMOVE B TRACKS FROM ITS FIT
          // ********************* todos los vertices primarios con constrain del Beam-Spot y escogemos el de mejor pointing angle ****************

        reco::Vertex vtxBSrf ;

        Double_t pVtxBSIPX_temp = -10000.0;
        Double_t pVtxBSIPY_temp = -10000.0;
        Double_t pVtxBSIPZ_temp = -10000.0;
        Double_t pVtxBSIPXE_temp = -10000.0;
        Double_t pVtxBSIPYE_temp = -10000.0;
        Double_t pVtxBSIPZE_temp = -10000.0;
        Double_t pVtxBSIPCL_temp = -10000.0;
        Double_t pVtxBSIPdN_temp = 0;
        Double_t lip = -100000.0;
        for(size_t i = 0; i < recVtxs->size(); ++i)
        {
            Double_t ptsum_ = 0;
            const Vertex &vtxBS = (*recVtxs)[i];
            vector<reco::TransientTrack> vertexTracks;
            for ( std::vector<TrackBaseRef >::const_iterator iTrack = vtxBS.tracks_begin(); iTrack != vtxBS.tracks_end(); ++iTrack)
            {
                TrackRef trackRef = iTrack->castTo<TrackRef>();
                // the  tracks in the D/B cand are
                if (  !( (iTrack1  ==trackRef.get() )
                       ||(vee->daughter(0)->bestTrack() == trackRef.get() )
                       ||(vee->daughter(1)->bestTrack() == trackRef.get() )) )
                {
                    // TransientTrack tt(trackRef, &(*bFieldHandle) );
                    TransientTrack tt = theB.build( fix_track(trackRef) );
                    vertexTracks.push_back(tt);
                    ptsum_ += trackRef->pt();
                } //else { std::cout << "found track match with primary" << endl;}
            }
            // if no tracks in primary or no reco track included in primary then don't do anything
            vtxBSrf = vtxBS;
            GlobalPoint PVRfP = GlobalPoint( vtxBS.x(), vtxBS.y(), vtxBS.z() );
            if ( vertexTracks.size()>0 && (vtxBS.tracksSize()!=vertexTracks.size()) )
            {
                AdaptiveVertexFitter theFitter;
                TransientVertex v = theFitter.vertex(vertexTracks,PVRfP);
                if ( v.isValid() )
                {
                    vtxBSrf = reco::Vertex(v);
                }
            }
            Double_t dx = (*Xi_vtx).position().x() - vtxBSrf.x();
            Double_t dy = (*Xi_vtx).position().y() - vtxBSrf.y();
            Double_t dz = (*Xi_vtx).position().z() - vtxBSrf.z();
            Double_t cosAlphaXYb = ( Xi_vFit->currentState().globalMomentum().x() * dx + Xi_vFit->currentState().globalMomentum().y()*dy + Xi_vFit->currentState().globalMomentum().z()*dz  )/( sqrt(dx*dx+dy*dy+dz*dz)* Xi_vFit->currentState().globalMomentum().mag() );
            if(cosAlphaXYb>lip)
            {
                lip = cosAlphaXYb ;
                pVtxBSIPX_temp     = vtxBSrf.x();
                pVtxBSIPY_temp     = vtxBSrf.y();
                pVtxBSIPZ_temp     = vtxBSrf.z();
                pVtxBSIPXE_temp    = vtxBSrf.covariance(0, 0);
                pVtxBSIPYE_temp    = vtxBSrf.covariance(1, 1);
                pVtxBSIPZE_temp    = vtxBSrf.covariance(2, 2);
                pVtxBSIPCL_temp    = (TMath::Prob(vtxBSrf.chi2(), (int)vtxBSrf.ndof()) );
                pVtxBSIPdN_temp    = vtxBS.tracksSize() - vertexTracks.size();
                bestVtxBSIP = vtxBSrf;
            }
        }
         // }}}

        // if ((Xi_mass_tmp > 1.28) && (Xi_mass_tmp < 1.29))
        // {
        // }
        //
        // topological requirements
        tmp_dx   = Xi_vtx->position().x() - pVtxBSIPX_temp ;
        tmp_dy   = Xi_vtx->position().y() - pVtxBSIPY_temp;
        tmp_edx  = fabs(Xi_vtx->error().cxx()) + fabs(pVtxBSIPXE_temp) ;
        tmp_edy  = fabs(Xi_vtx->error().cyy()) + fabs(pVtxBSIPYE_temp) ;
        tmp_px   = Xi_vFit->currentState().globalMomentum().x();
        tmp_py   = Xi_vFit->currentState().globalMomentum().y();
        //
        tmp_DS2  = sqrt((tmp_dx*tmp_dx)/tmp_edx + (tmp_dy*tmp_dy)/tmp_edy );
        if( tmp_DS2 < 2) continue;
        tmp_Cos2 = (tmp_dx * tmp_px + tmp_dy * tmp_py) / (sqrt(tmp_dx * tmp_dx + tmp_dy * tmp_dy)*sqrt(tmp_px * tmp_px + tmp_py * tmp_py ));
        if( tmp_Cos2 < 0.5) continue;
        // //
        //
        // if (fabs(iTrack1->dxy(bestVtxBSIP.position())) / (0.000001 + fabs(iTrack1->dxyError()) ) < 0.5 ) continue;
        //
        // extra checks that tracks arent the same
        if (( fabs(vee->daughter(0)->bestTrack()->dxy(bestVtxBSIP.position())) / (0.000001 + fabs(vee->daughter(0)->bestTrack()->dxyError())) == fabs(iTrack1->dxy(bestVtxBSIP.position())) / (0.000001 + fabs(iTrack1->dxyError())) ) && ( iTrack1->numberOfValidHits() == vee->daughter(0)->bestTrack()->numberOfValidHits() )) continue;
        if (( fabs(vee->daughter(1)->bestTrack()->dxy(bestVtxBSIP.position())) / (0.000001 + fabs(vee->daughter(1)->bestTrack()->dxyError())) == fabs(iTrack1->dxy(bestVtxBSIP.position())) / (0.000001 + fabs(iTrack1->dxyError())) ) && ( iTrack1->numberOfValidHits() == vee->daughter(1)->bestTrack()->numberOfValidHits() )) continue;
        //
        LA_mass         ->push_back( MY_lamass2                                             );
        LA_vtxp         ->push_back( LA_Prob_tmp                                            );
        LA_pr_charg     ->push_back( protTT.charge()                                        );
        LA_px_CL        ->push_back( vee_vFit_withMC->currentState().globalMomentum().x()   );
        LA_py_CL        ->push_back( vee_vFit_withMC->currentState().globalMomentum().y()   );
        LA_pz_CL        ->push_back( vee_vFit_withMC->currentState().globalMomentum().z()   );
        LA_VtxX_CL      ->push_back( veeVtx_vFit_withMC->position().x()                     );
        LA_VtxY_CL      ->push_back( veeVtx_vFit_withMC->position().y()                     );
        LA_VtxZ_CL      ->push_back( veeVtx_vFit_withMC->position().z()                     );
        LA_VtxXE_CL     ->push_back( veeVtx_vFit_withMC->error().cxx()                      );
        LA_VtxYE_CL     ->push_back( veeVtx_vFit_withMC->error().cyy()                      );
        LA_VtxZE_CL     ->push_back( veeVtx_vFit_withMC->error().czz()                      );
        LA_pi_px        ->push_back( pionCand->currentState().globalMomentum().x()          );
        LA_pi_py        ->push_back( pionCand->currentState().globalMomentum().y()          );
        LA_pi_pz        ->push_back( pionCand->currentState().globalMomentum().z()          );
        LA_pr_px        ->push_back( protCand->currentState().globalMomentum().x()          );
        LA_pr_py        ->push_back( protCand->currentState().globalMomentum().y()          );
        LA_pr_pz        ->push_back( protCand->currentState().globalMomentum().z()          );
        LA_pi_ips       ->push_back( fabs(vee->daughter(0)->bestTrack()->dxy(bestVtxBSIP.position())) / (0.000001 + fabs(vee->daughter(0)->bestTrack()->dxyError())) );
        //
        XI_mass         ->push_back( Xi_mass_tmp                                            );
        XI_vtxp         ->push_back( Xi_Prob_tmp                                            );
        XI_pi_charg     ->push_back( iTrack1->charge()                                      );
        XI_pi_purit     ->push_back( pi1quality                                             );
        XI_pi_lost      ->push_back( 0                                                      );
        XI_px           ->push_back( Xi_vFit->currentState().globalMomentum().x()           );
        XI_py           ->push_back( Xi_vFit->currentState().globalMomentum().y()           );
        XI_pz           ->push_back( Xi_vFit->currentState().globalMomentum().z()           );
        XI_pi_px        ->push_back( p4pi1.Px()                                             );
        XI_pi_py        ->push_back( p4pi1.Py()                                             );
        XI_pi_pz        ->push_back( p4pi1.Pz()                                             );
        XI_pi_px_CV     ->push_back( pixiCand->currentState().globalMomentum().x()          );
        XI_pi_py_CV     ->push_back( pixiCand->currentState().globalMomentum().y()          );
        XI_pi_pz_CV     ->push_back( pixiCand->currentState().globalMomentum().z()          );
        XI_VtxX         ->push_back( Xi_vtx->position().x()                                 );
        XI_VtxY         ->push_back( Xi_vtx->position().y()                                 );
        XI_VtxZ         ->push_back( Xi_vtx->position().z()                                 );
        XI_VtxXE        ->push_back( Xi_vtx->error().cxx()                                  );
        XI_VtxYE        ->push_back( Xi_vtx->error().cyy()                                  );
        XI_VtxZE        ->push_back( Xi_vtx->error().czz()                                  );
        XI_pi_hit       ->push_back( iTrack1->numberOfValidHits()                 );
        XI_pi_pix       ->push_back( iTrack1->hitPattern().numberOfValidPixelHits());
        XI_pi_ips       ->push_back( fabs(iTrack1->dxy(bestVtxBSIP.position())) / (0.000001 + fabs(iTrack1->dxyError())) );
        ////
        PV_becos_XX     ->push_back( pVtxBSIPX_temp                                         );
        PV_becos_YY     ->push_back( pVtxBSIPY_temp                                         );
        PV_becos_ZZ     ->push_back( pVtxBSIPZ_temp                                         );
        PV_becos_EX     ->push_back( pVtxBSIPXE_temp                                        );
        PV_becos_EY     ->push_back( pVtxBSIPYE_temp                                        );
        PV_becos_EZ     ->push_back( pVtxBSIPZE_temp                                        );
        PV_becos_CL     ->push_back( pVtxBSIPCL_temp                                        );
        PV_becos_dN     ->push_back( pVtxBSIPdN_temp                                        );
        //
        nCand ++;
        Xi_candidate_init   .clear();
    }
    //
    if (numLost < 1 ) continue;
    int lotnumb  = 0;
    /// LOOPING OVER LOST TRACKS IN MINIAOD
    for (pat::PackedCandidateCollection::const_iterator iTra1 = lostTracks->begin(); iTra1 != lostTracks->end(); iTra1++)
    {
        // std::cout << " new lost trk sss.. "<< std::endl;
        if (!(iTra1->hasTrackDetails())) continue;
        auto iTrack1 = iTra1->bestTrack();
        if (iTrack1 == nullptr) continue;
        //
        // std::cout << " new lost trk 0.. " << run << "," << event <<  " " << iTrack1->pt() << std::endl;
        //
        lotids.insert(lotnumb);
        lotnumb += 1;
        //
        // if(iTrack1->pt()<0.2) continue;
        if( fabs(iTrack1->eta()) > 3.0 ) continue;
        //
        if(iTrack1 == vee->daughter(0)->bestTrack() ) continue;
        if(iTrack1 == vee->daughter(0)->bestTrack() ) continue;
        if (iTrack1->charge() * protTT.charge() > 0) continue;
        // if (!(iTrack1->quality(reco::TrackBase::highPurity))) continue;
        int pi1quality = 0; // 1 = loose, 2 = tight, 4 = highPurity
        if(iTrack1->quality(reco::TrackBase::loose))
        {
            pi1quality += 1;
        }
        if(iTrack1->quality(reco::TrackBase::tight))
        {
            pi1quality += 2;
        }
        if(iTrack1->quality(reco::TrackBase::highPurity))
        {
            pi1quality += 4;
        }
        //
        fitgood = 1;
        TransientTrack pion1TT;
        try
        {
            TransientTrack pion1TT = theB.build( fix_track(iTrack1) );
            if (!pion1TT.isValid()) fitgood = 0;
        }
        catch (const Exception &)
        {
            std::cout << " CRASH building TT at run evt " << run << "," << event << std::endl;
            fitgood = 0;
        }
        if (fitgood == 0) continue;
        // std::cout << " before valid check 1.. " << run << "," << event <<  " " << iTrack1->pt() << std::endl;
        // if (!pion1TT.isValid()) continue;
        // std::cout << " after valid check 2.. " << run << "," << event <<  " " << std::endl;
        fitgood = 1;
        try
        {
            if (pionTT == pion1TT) fitgood = 0;
        }
        catch (const Exception &)
        {
            std::cout << " CRASH checking TT 2.5 at run evt " << run << "," << event << std::endl;
            fitgood = 0;
        }
        if (fitgood == 0) continue;
        //
        // std::cout << " after check 3.. " << run << "," << event <<  " " << std::endl;
        //
        TLorentzVector p4pi1;
        p4pi1.SetPtEtaPhiM(iTrack1->pt(),iTrack1->eta(),iTrack1->phi(), PDG_PION_MASS);
        //
        // HERE COMPUTING THE STUFF
        if ((p4pi1 + p4la).M() < 1.45)
        {
            loXids.insert(lotnumb);
        }
        //
        //
        if ((p4pi1 + p4la).M() > 2.0) continue; //
        //
        // std::cout << " before valid check 5.. " << run << "," << event <<  " " << std::endl;
        fitgood = 1;
        try
        {
            if (!pion1TT.isValid()) fitgood = 0;
        }
        catch (const Exception &)
        {
            std::cout << " CRASH checking valid 6 TT at run evt " << run << "," << event << std::endl;
            fitgood = 0;
        }
        if (fitgood == 0) continue;
        //
        std::vector<RefCountedKinematicParticle> Xi_candidate_init;
        fitgood = 1;
        try
        {
            Xi_candidate_init.push_back(pFactory.particle(pion1TT, PM_PDG_PION_MASS, chi,ndf, PM_sigma));
        }
        catch (const Exception &)
        {
            std::cout << " CRASH at run evt 6p5 " << run << "," << event << std::endl;
            fitgood = 0;
        }
        if (fitgood == 0) continue;
        Xi_candidate_init.push_back(vee_vFit_withMC);
        RefCountedKinematicTree XiTree;
        KinematicParticleVertexFitter xFitter; //KinematicParticleVertexFitter
        fitgood = 1;
        try
        {
            XiTree = xFitter.fit(Xi_candidate_init);
        }
        catch (const VertexException &)
        {
            fitgood = 0;
        }
        if (fitgood == 0) continue;
        //
        if (!XiTree->isValid()) continue;
        XiTree->movePointerToTheTop();
        RefCountedKinematicParticle Xi_vFit = XiTree->currentParticle();
        RefCountedKinematicVertex   Xi_vtx = XiTree->currentDecayVertex();
        //
        double Xi_mass_tmp = Xi_vFit->currentState().mass();
        if (Xi_mass_tmp > 1.36) continue; //
        if (Xi_vtx->chiSquared() < 0) continue;
        Xi_Prob_tmp   = (double) (0.0 + TMath::Prob(Xi_vtx->chiSquared(), Xi_vtx->degreesOfFreedom()));
        if(Xi_Prob_tmp < 0.01) continue;
        //
        XiTree->movePointerToTheFirstChild();
        RefCountedKinematicParticle pixiCand    = XiTree->currentParticle();
        XiTree->movePointerToTheTop();
        //
        // std::cout << " Good VTX Fit with Lost Trk at run evt " << run << "," << event << std::endl;
        //
        // topological requirements
        double tmp_dx   = veeVtx_vFit_withMC->position().x() - Xi_vtx->position().x() ;
        double tmp_dy   = veeVtx_vFit_withMC->position().y() - Xi_vtx->position().y() ;
        double tmp_edx  = fabs(veeVtx_vFit_withMC->error().cxx()) + fabs(Xi_vtx->error().cxx()) ;
        double tmp_edy  = fabs(veeVtx_vFit_withMC->error().cyy()) + fabs(Xi_vtx->error().cyy()) ;
        double tmp_px   = vee_vFit_withMC->currentState().globalMomentum().x();
        double tmp_py   = vee_vFit_withMC->currentState().globalMomentum().y();
        //
        double tmp_DS2  = sqrt((tmp_dx*tmp_dx)/tmp_edx + (tmp_dy*tmp_dy)/tmp_edy );
        if( tmp_DS2 < 2) continue;
        double tmp_Cos2 = (tmp_dx * tmp_px + tmp_dy * tmp_py) / (sqrt(tmp_dx * tmp_dx + tmp_dy * tmp_dy)*sqrt(tmp_px * tmp_px + tmp_py * tmp_py ));
        if( tmp_Cos2 < 0.5) continue;
        //
        ////
        reco::Vertex bestVtxBSIP;  /// using Xi_vFit and Xi_vtx HERE !!
         // {{{ GET THE BEST PV BY CHOSING THE BEST POINTING ANGLE AND REMOVE B TRACKS FROM ITS FIT
          // ********************* todos los vertices primarios con constrain del Beam-Spot y escogemos el de mejor pointing angle ****************

        reco::Vertex vtxBSrf ;

        Double_t pVtxBSIPX_temp = -10000.0;
        Double_t pVtxBSIPY_temp = -10000.0;
        Double_t pVtxBSIPZ_temp = -10000.0;
        Double_t pVtxBSIPXE_temp = -10000.0;
        Double_t pVtxBSIPYE_temp = -10000.0;
        Double_t pVtxBSIPZE_temp = -10000.0;
        Double_t pVtxBSIPCL_temp = -10000.0;
        Double_t pVtxBSIPdN_temp = 0;
        Double_t lip = -100000.0;
        for(size_t i = 0; i < recVtxs->size(); ++i)
        {
            Double_t ptsum_ = 0;
            const Vertex &vtxBS = (*recVtxs)[i];
            vector<reco::TransientTrack> vertexTracks;
            for ( std::vector<TrackBaseRef >::const_iterator iTrack = vtxBS.tracks_begin(); iTrack != vtxBS.tracks_end(); ++iTrack)
            {
                TrackRef trackRef = iTrack->castTo<TrackRef>();
                // the  tracks in the D/B cand are
                if (  !( (iTrack1  ==trackRef.get() )
                       ||(vee->daughter(0)->bestTrack() == trackRef.get() )
                       ||(vee->daughter(1)->bestTrack() == trackRef.get() )) )
                {
                    // TransientTrack tt(trackRef, &(*bFieldHandle) );
                    TransientTrack tt = theB.build( fix_track(trackRef) );
                    vertexTracks.push_back(tt);
                    ptsum_ += trackRef->pt();
                } //else { std::cout << "found track match with primary" << endl;}
            }
            // if no tracks in primary or no reco track included in primary then don't do anything
            vtxBSrf = vtxBS;
            GlobalPoint PVRfP = GlobalPoint( vtxBS.x(), vtxBS.y(), vtxBS.z() );
            if ( vertexTracks.size()>0 && (vtxBS.tracksSize()!=vertexTracks.size()) )
            {
                AdaptiveVertexFitter theFitter;
                TransientVertex v = theFitter.vertex(vertexTracks,PVRfP);
                if ( v.isValid() )
                {
                    vtxBSrf = reco::Vertex(v);
                }
            }
            Double_t dx = (*Xi_vtx).position().x() - vtxBSrf.x();
            Double_t dy = (*Xi_vtx).position().y() - vtxBSrf.y();
            Double_t dz = (*Xi_vtx).position().z() - vtxBSrf.z();
            Double_t cosAlphaXYb = ( Xi_vFit->currentState().globalMomentum().x() * dx + Xi_vFit->currentState().globalMomentum().y()*dy + Xi_vFit->currentState().globalMomentum().z()*dz  )/( sqrt(dx*dx+dy*dy+dz*dz)* Xi_vFit->currentState().globalMomentum().mag() );
            if(cosAlphaXYb>lip)
            {
                lip = cosAlphaXYb ;
                pVtxBSIPX_temp     = vtxBSrf.x();
                pVtxBSIPY_temp     = vtxBSrf.y();
                pVtxBSIPZ_temp     = vtxBSrf.z();
                pVtxBSIPXE_temp    = vtxBSrf.covariance(0, 0);
                pVtxBSIPYE_temp    = vtxBSrf.covariance(1, 1);
                pVtxBSIPZE_temp    = vtxBSrf.covariance(2, 2);
                pVtxBSIPCL_temp    = (TMath::Prob(vtxBSrf.chi2(), (int)vtxBSrf.ndof()) );
                pVtxBSIPdN_temp    = vtxBS.tracksSize() - vertexTracks.size();
                bestVtxBSIP = vtxBSrf;
            }
        }
         // }}}

        // if ((Xi_mass_tmp > 1.28) && (Xi_mass_tmp < 1.29))
        // {
        // }
        //
        // topological requirements
        tmp_dx   = Xi_vtx->position().x() - pVtxBSIPX_temp ;
        tmp_dy   = Xi_vtx->position().y() - pVtxBSIPY_temp;
        tmp_edx  = fabs(Xi_vtx->error().cxx()) + fabs(pVtxBSIPXE_temp) ;
        tmp_edy  = fabs(Xi_vtx->error().cyy()) + fabs(pVtxBSIPYE_temp) ;
        tmp_px   = Xi_vFit->currentState().globalMomentum().x();
        tmp_py   = Xi_vFit->currentState().globalMomentum().y();
        //
        tmp_DS2  = sqrt((tmp_dx*tmp_dx)/tmp_edx + (tmp_dy*tmp_dy)/tmp_edy );
        if( tmp_DS2 < 2) continue;
        tmp_Cos2 = (tmp_dx * tmp_px + tmp_dy * tmp_py) / (sqrt(tmp_dx * tmp_dx + tmp_dy * tmp_dy)*sqrt(tmp_px * tmp_px + tmp_py * tmp_py ));
        if( tmp_Cos2 < 0.5) continue;
        // //
        //
        // if (fabs(iTrack1->dxy(bestVtxBSIP.position())) / (0.000001 + fabs(iTrack1->dxyError()) ) < 0.5 ) continue;
        //
        // extra checks that tracks arent the same
        if (( fabs(vee->daughter(0)->bestTrack()->dxy(bestVtxBSIP.position())) / (0.000001 + fabs(vee->daughter(0)->bestTrack()->dxyError())) == fabs(iTrack1->dxy(bestVtxBSIP.position())) / (0.000001 + fabs(iTrack1->dxyError())) ) && ( iTrack1->numberOfValidHits() == vee->daughter(0)->bestTrack()->numberOfValidHits() )) continue;
        if (( fabs(vee->daughter(1)->bestTrack()->dxy(bestVtxBSIP.position())) / (0.000001 + fabs(vee->daughter(1)->bestTrack()->dxyError())) == fabs(iTrack1->dxy(bestVtxBSIP.position())) / (0.000001 + fabs(iTrack1->dxyError())) ) && ( iTrack1->numberOfValidHits() == vee->daughter(1)->bestTrack()->numberOfValidHits() )) continue;
        //
        LA_mass         ->push_back( MY_lamass2                                             );
        LA_vtxp         ->push_back( LA_Prob_tmp                                            );
        LA_pr_charg     ->push_back( protTT.charge()                                        );
        LA_px_CL        ->push_back( vee_vFit_withMC->currentState().globalMomentum().x()   );
        LA_py_CL        ->push_back( vee_vFit_withMC->currentState().globalMomentum().y()   );
        LA_pz_CL        ->push_back( vee_vFit_withMC->currentState().globalMomentum().z()   );
        LA_VtxX_CL      ->push_back( veeVtx_vFit_withMC->position().x()                     );
        LA_VtxY_CL      ->push_back( veeVtx_vFit_withMC->position().y()                     );
        LA_VtxZ_CL      ->push_back( veeVtx_vFit_withMC->position().z()                     );
        LA_VtxXE_CL     ->push_back( veeVtx_vFit_withMC->error().cxx()                      );
        LA_VtxYE_CL     ->push_back( veeVtx_vFit_withMC->error().cyy()                      );
        LA_VtxZE_CL     ->push_back( veeVtx_vFit_withMC->error().czz()                      );
        LA_pi_px        ->push_back( pionCand->currentState().globalMomentum().x()          );
        LA_pi_py        ->push_back( pionCand->currentState().globalMomentum().y()          );
        LA_pi_pz        ->push_back( pionCand->currentState().globalMomentum().z()          );
        LA_pr_px        ->push_back( protCand->currentState().globalMomentum().x()          );
        LA_pr_py        ->push_back( protCand->currentState().globalMomentum().y()          );
        LA_pr_pz        ->push_back( protCand->currentState().globalMomentum().z()          );
        LA_pi_ips       ->push_back( fabs(vee->daughter(0)->bestTrack()->dxy(bestVtxBSIP.position())) / (0.000001 + fabs(vee->daughter(0)->bestTrack()->dxyError())) );
        //
        XI_mass         ->push_back( Xi_mass_tmp                                            );
        XI_vtxp         ->push_back( Xi_Prob_tmp                                            );
        XI_pi_charg     ->push_back( iTrack1->charge()                                      );
        XI_pi_purit     ->push_back( pi1quality                                             );
        XI_pi_lost      ->push_back( 1                                                      );
        XI_px           ->push_back( Xi_vFit->currentState().globalMomentum().x()           );
        XI_py           ->push_back( Xi_vFit->currentState().globalMomentum().y()           );
        XI_pz           ->push_back( Xi_vFit->currentState().globalMomentum().z()           );
        XI_pi_px        ->push_back( p4pi1.Px()                                             );
        XI_pi_py        ->push_back( p4pi1.Py()                                             );
        XI_pi_pz        ->push_back( p4pi1.Pz()                                             );
        XI_pi_px_CV     ->push_back( pixiCand->currentState().globalMomentum().x()          );
        XI_pi_py_CV     ->push_back( pixiCand->currentState().globalMomentum().y()          );
        XI_pi_pz_CV     ->push_back( pixiCand->currentState().globalMomentum().z()          );
        XI_VtxX         ->push_back( Xi_vtx->position().x()                                 );
        XI_VtxY         ->push_back( Xi_vtx->position().y()                                 );
        XI_VtxZ         ->push_back( Xi_vtx->position().z()                                 );
        XI_VtxXE        ->push_back( Xi_vtx->error().cxx()                                  );
        XI_VtxYE        ->push_back( Xi_vtx->error().cyy()                                  );
        XI_VtxZE        ->push_back( Xi_vtx->error().czz()                                  );
        XI_pi_hit       ->push_back( iTrack1->numberOfValidHits()                 );
        XI_pi_pix       ->push_back( iTrack1->hitPattern().numberOfValidPixelHits());
        XI_pi_ips       ->push_back( fabs(iTrack1->dxy(bestVtxBSIP.position())) / (0.000001 + fabs(iTrack1->dxyError())) );
        ////
        PV_becos_XX     ->push_back( pVtxBSIPX_temp                                         );
        PV_becos_YY     ->push_back( pVtxBSIPY_temp                                         );
        PV_becos_ZZ     ->push_back( pVtxBSIPZ_temp                                         );
        PV_becos_EX     ->push_back( pVtxBSIPXE_temp                                        );
        PV_becos_EY     ->push_back( pVtxBSIPYE_temp                                        );
        PV_becos_EZ     ->push_back( pVtxBSIPZE_temp                                        );
        PV_becos_CL     ->push_back( pVtxBSIPCL_temp                                        );
        PV_becos_dN     ->push_back( pVtxBSIPdN_temp                                        );
        //
        nCand ++;
        Xi_candidate_init   .clear();
    }
    //
    VeeParticles        .clear();
}

// numTrack = trkids.size();
numLosa = lotids.size();
numLosx = loXids.size();

// for (reco::MuonCollection::const_iterator muon=muonColl->begin(), muonCollEnd=muonColl->end(); muon!=muonCollEnd;  ++muon)
// {
//     if(!(
//     muon->isGlobalMuon() &&
//     muon->globalTrack()->normalizedChi2() < 10 &&
//     muon->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 &&
//     muon->numberOfMatchedStations() > 1 &&
//     muon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0 &&
//     muon::isGoodMuon(*muon, muon::TMLastStationTight) &&
//     muon->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 4 &&
//     muon->pt() > 3.0 && fabs(muon->eta()) < 2.4
//     )) continue;
//
//     mu_pt->push_back( muon->pt() ) ;
//     nCand ++;
// }

// ===================== END OF EVENT : WRITE ETC ++++++++++++++++++++++

// if (nCand > 0)
// {
    wwtree->Fill();
// } // nCand > 0

//
LA_mass->clear();       LA_vtxp->clear();       LA_pr_charg->clear();
LA_px_CL->clear();      LA_py_CL->clear();      LA_pz_CL->clear();
LA_VtxX_CL ->clear();   LA_VtxY_CL ->clear();   LA_VtxZ_CL ->clear();
LA_VtxXE_CL->clear();   LA_VtxYE_CL->clear();   LA_VtxZE_CL->clear();
LA_pi_px->clear();      LA_pi_py->clear();      LA_pi_pz->clear();
LA_pr_px->clear();      LA_pr_py->clear();      LA_pr_pz->clear();
LA_pi_ips->clear();
//
XI_mass->clear();       XI_vtxp->clear();       XI_pi_charg->clear();
XI_pi_purit->clear();   XI_pi_lost->clear();
XI_px->clear();         XI_py->clear();         XI_pz->clear();
XI_pi_px->clear();      XI_pi_py->clear();      XI_pi_pz->clear();
XI_pi_px_CV->clear();   XI_pi_py_CV->clear();   XI_pi_pz_CV->clear();
XI_VtxX ->clear();      XI_VtxY ->clear();      XI_VtxZ ->clear();
XI_VtxXE->clear();      XI_VtxYE->clear();      XI_VtxZE->clear();
XI_pi_hit->clear();     XI_pi_pix->clear();     XI_pi_ips->clear();

PV_becos_XX->clear();   PV_becos_YY->clear();   PV_becos_ZZ->clear();
PV_becos_EX->clear();   PV_becos_EY->clear();   PV_becos_EZ->clear();
PV_becos_CL->clear();   PV_becos_dN->clear();

gen_XI_px->clear();     gen_XI_py->clear();     gen_XI_pz->clear(); gen_XI_ma->clear();
}


// ------------ method called once each job just before starting event loop  ------------
void XiAnalyzer::beginJob()
{
using namespace std;
using namespace reco;
//
cout << "------------------------------->>>>> Begin Job" << endl;

f = new TFile(fileName.c_str(), "RECREATE");
wwtree  = new TTree("wztree", "muons tree");

wwtree->Branch("nCand"              , &nCand            , "nCand/I"     );

wwtree->Branch("run"                , &run              , "run/I"       );
wwtree->Branch("event"              , &event            , "event/I"     );
wwtree->Branch("lumi"               , &lumi             , "lumi/F"      );

wwtree->Branch("numPV"              , &numPV            , "numPV/I"     );
wwtree->Branch("numTrack"           , &numTrack         , "numTrack/I"  );
wwtree->Branch("numLost"            , &numLost          , "numLost/I"   );
wwtree->Branch("numLosa"            , &numLosa          , "numLosa/I"   );
wwtree->Branch("numLosx"            , &numLosx          , "numLosx/I"   );
wwtree->Branch("numV0"              , &numV0            , "numV0/I"     );

wwtree->Branch("LA_mass"            , &LA_mass          );
wwtree->Branch("LA_vtxp"            , &LA_vtxp          );
wwtree->Branch("LA_pr_charg"        , &LA_pr_charg      );
wwtree->Branch("LA_px_CL"           , &LA_px_CL         );
wwtree->Branch("LA_py_CL"           , &LA_py_CL         );
wwtree->Branch("LA_pz_CL"           , &LA_pz_CL         );
wwtree->Branch("LA_VtxX_CL"         , &LA_VtxX_CL       );
wwtree->Branch("LA_VtxY_CL"         , &LA_VtxY_CL       );
wwtree->Branch("LA_VtxZ_CL"         , &LA_VtxZ_CL       );
wwtree->Branch("LA_VtxXE_CL"        , &LA_VtxXE_CL      );
wwtree->Branch("LA_VtxYE_CL"        , &LA_VtxYE_CL      );
wwtree->Branch("LA_VtxZE_CL"        , &LA_VtxZE_CL      );
wwtree->Branch("LA_pi_px"           , &LA_pi_px         );
wwtree->Branch("LA_pi_py"           , &LA_pi_py         );
wwtree->Branch("LA_pi_pz"           , &LA_pi_pz         );
wwtree->Branch("LA_pr_px"           , &LA_pr_px         );
wwtree->Branch("LA_pr_py"           , &LA_pr_py         );
wwtree->Branch("LA_pr_pz"           , &LA_pr_pz         );
wwtree->Branch("LA_pi_ips"          , &LA_pi_ips        );

wwtree->Branch("XI_mass"            , &XI_mass          );
wwtree->Branch("XI_vtxp"            , &XI_vtxp          );
wwtree->Branch("XI_pi_charg"        , &XI_pi_charg      );
wwtree->Branch("XI_pi_purit"        , &XI_pi_purit      );
wwtree->Branch("XI_pi_lost"         , &XI_pi_lost       );
wwtree->Branch("XI_px"              , &XI_px            );
wwtree->Branch("XI_py"              , &XI_py            );
wwtree->Branch("XI_pz"              , &XI_pz            );
wwtree->Branch("XI_pi_px"           , &XI_pi_px         );
wwtree->Branch("XI_pi_py"           , &XI_pi_py         );
wwtree->Branch("XI_pi_pz"           , &XI_pi_pz         );
wwtree->Branch("XI_pi_px_CV"        , &XI_pi_px_CV      );
wwtree->Branch("XI_pi_py_CV"        , &XI_pi_py_CV      );
wwtree->Branch("XI_pi_pz_CV"        , &XI_pi_pz_CV      );
wwtree->Branch("XI_VtxX"            , &XI_VtxX          );
wwtree->Branch("XI_VtxY"            , &XI_VtxY          );
wwtree->Branch("XI_VtxZ"            , &XI_VtxZ          );
wwtree->Branch("XI_VtxXE"           , &XI_VtxXE         );
wwtree->Branch("XI_VtxYE"           , &XI_VtxYE         );
wwtree->Branch("XI_VtxZE"           , &XI_VtxZE         );
wwtree->Branch("XI_pi_hit"          , &XI_pi_hit        );
wwtree->Branch("XI_pi_pix"          , &XI_pi_pix        );
wwtree->Branch("XI_pi_ips"          , &XI_pi_ips        );

wwtree->Branch("PV_becos_XX"        , &PV_becos_XX      );
wwtree->Branch("PV_becos_YY"        , &PV_becos_YY      );
wwtree->Branch("PV_becos_ZZ"        , &PV_becos_ZZ      );
wwtree->Branch("PV_becos_EX"        , &PV_becos_EX      );
wwtree->Branch("PV_becos_EY"        , &PV_becos_EY      );
wwtree->Branch("PV_becos_EZ"        , &PV_becos_EZ      );
wwtree->Branch("PV_becos_CL"        , &PV_becos_CL      );
wwtree->Branch("PV_becos_dN"        , &PV_becos_dN      );

wwtree->Branch("gen_XI_px"          , &gen_XI_px        );
wwtree->Branch("gen_XI_py"          , &gen_XI_py        );
wwtree->Branch("gen_XI_pz"          , &gen_XI_pz        );
wwtree->Branch("gen_XI_ma"          , &gen_XI_ma        );
// wwtree1  = new TTree("gen", "muons tree");
// wwtree1->Branch("gen_xim_p4",   "TLorentzVector",   &gen_xim_p4);

}

// ------------ method called once each job just after ending the event loop  ------------
void
XiAnalyzer::endJob()
{

using namespace std;
cout << "------------------------------->>>>> End Job" << endl;
f->WriteTObject(wwtree);
delete wwtree;
f->Close();

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
XiAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(XiAnalyzer);
