// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>

#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TRandom.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

//
// class decleration
//

#define PI 3.1416

using namespace std;

class D0Corr : public edm::EDAnalyzer {
public:
    explicit D0Corr(const edm::ParameterSet&);
    ~D0Corr();
    
    
private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;
    
    // ----------member data ---------------------------
    
    TH1D* hMult;
    
    TH2F* effhisto;
    TH2D* effhisto_d0;
    
    TH1D* hMass_d0[13];
    TH1D* hKET_d0[13];
    TH1D* hPt_d0[13];
    TH1D* hEta_d0[13];
    TH1D* hMult_d0[13];
    TH2D* hSignal_d0[13];
    TH2D* hBackground_d0[13];
    
    TH1D* hKET_d0_bkg[13];
    TH1D* hPt_d0_bkg[13];
    TH1D* hEta_d0_bkg[13];
    TH1D* hMult_d0_bkg[13];
    TH2D* hSignal_d0_bkg[13];
    TH2D* hBackground_d0_bkg[13];
    
    TH1D* hMult_ass;
    
    vector<TVector3> *pVect_trg_d0[13];
    vector< vector<TVector3> > *pVectVect_trg_d0[13];
    
    vector<TVector3> *pVect_dau_d0[13];
    vector< vector<TVector3> > *pVectVect_dau_d0[13];
    
    vector<TVector3> *pVect_trg_d0_bkg[13];
    vector< vector<TVector3> > *pVectVect_trg_d0_bkg[13];
    
    vector<TVector3> *pVect_dau_d0_bkg[13];
    vector< vector<TVector3> > *pVectVect_dau_d0_bkg[13];
    
    vector<TVector3> *pVect_ass;
    vector< vector<TVector3> > *pVectVect_ass;
    vector<double> *zvtxVect;
    
    vector<double> ptcut_d0_;
    vector<double> sigma_d0_;
    vector<double> mean_d0_;
    vector<double> angle_d0_;
    vector<double> dls_d0_;
    vector<double> vtxprob_d0_;
    vector<double> daupt_d0_;
    
    int ptbin_n_;
    double etaMin_trg_;
    double etaMax_trg_;
    double etaMin_ass_;
    double etaMax_ass_;
    double ptMin_ass_;
    double ptMax_ass_;
    double multMax_;
    double multMin_;
    int bkgFactor_;
    double peakFactor_;
    double sideFactor_;
    double mis_d0_range_;
    double mis_la_range_;
    double mis_ph_range_;
    bool rejectDaughter_;
    double d0_eta_;
    double d0_rapidity_;
    double d0_dau_eta_;
    int d0_dau_nhit_;
    double d0_dau_pterr_;

    
    edm::EDGetTokenT<reco::VertexCollection> tok_offlinePV_;
    edm::EDGetTokenT<reco::TrackCollection> tok_generalTrk_;
    
    edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> tok_d0_;
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

D0Corr::D0Corr(const edm::ParameterSet& iConfig)
{
    
    //now do what ever initialization is needed
    etaMin_trg_ = iConfig.getUntrackedParameter<double>("etaMin_trg", -2.4);
    etaMax_trg_ = iConfig.getUntrackedParameter<double>("etaMax_trg", 2.4);
    etaMin_ass_ = iConfig.getUntrackedParameter<double>("etaMin_ass", -2.4);
    etaMax_ass_ = iConfig.getUntrackedParameter<double>("etaMax_ass", 2.4);
    ptMin_ass_ = iConfig.getUntrackedParameter<double>("ptMin_ass", 0.3);
    ptMax_ass_ = iConfig.getUntrackedParameter<double>("ptMax_ass", 3.0);
    multMax_ = iConfig.getUntrackedParameter<double>("multMax", 220);
    multMin_ = iConfig.getUntrackedParameter<double>("multMin", 185);
    bkgFactor_ = iConfig.getUntrackedParameter<int>("bkgFactor", 20);
    ptbin_n_ = iConfig.getUntrackedParameter<int>("ptbin_n", 13);
    peakFactor_ = iConfig.getUntrackedParameter<double>("peakFactor", 2.0);
    sideFactor_ = iConfig.getUntrackedParameter<double>("sideFactor", 3.0);
    mis_d0_range_ = iConfig.getUntrackedParameter<double>("mis_d0_range", 0.020);
    mis_la_range_ = iConfig.getUntrackedParameter<double>("mis_la_range", 0.010);
    mis_ph_range_ = iConfig.getUntrackedParameter<double>("mis_ph_range", 0.015);
    d0_eta_ = iConfig.getUntrackedParameter<double>("d0_eta", 9999.0);
    d0_rapidity_ = iConfig.getUntrackedParameter<double>("d0_rapidity", 1.0);
    d0_dau_eta_ = iConfig.getUntrackedParameter<double>("d0_dau_eta", 1.5);
    d0_dau_nhit_ = iConfig.getUntrackedParameter<int>("d0_dau_nhit", 11);
    d0_dau_pterr_ = iConfig.getUntrackedParameter<double>("d0_dau_pterr", 0.1);

    
    sigma_d0_ = iConfig.getUntrackedParameter<std::vector<double> >("sigma_d0");
    mean_d0_ = iConfig.getUntrackedParameter<std::vector<double> >("mean_d0");
    ptcut_d0_ = iConfig.getUntrackedParameter<std::vector<double> >("ptcut_d0");
    
    angle_d0_ = iConfig.getUntrackedParameter<std::vector<double> >("angle_d0");
    dls_d0_ = iConfig.getUntrackedParameter<std::vector<double> >("dls_d0");
    vtxprob_d0_ = iConfig.getUntrackedParameter<std::vector<double> >("vtxprob_d0");
    daupt_d0_ = iConfig.getUntrackedParameter<std::vector<double> >("daupt_d0");


    rejectDaughter_ = iConfig.getUntrackedParameter<bool>("rejectDaughter");
    
    tok_offlinePV_ = consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"));
    tok_generalTrk_ = consumes<reco::TrackCollection>(edm::InputTag("generalTracks"));
    tok_d0_ = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("D0Collection"));
}


D0Corr::~D0Corr()
{
    
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
D0Corr::analyze(const edm::Event& iEvent, const edm::EventSetup&
                         iSetup)
{
    using namespace edm;
    
    // select on requirement of valid vertex
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(tok_offlinePV_,vertices);
    double bestvz=-999.9, bestvx=-999.9, bestvy=-999.9;
    double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
    const reco::Vertex & vtx = (*vertices)[0];
    bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
    bestvzError = vtx.zError(); bestvxError = vtx.xError(); bestvyError = vtx.yError();
    
    if(bestvz < -15.0 || bestvz>15.0) return;

    edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates_d0;
    iEvent.getByToken(tok_d0_,v0candidates_d0);
    if(!v0candidates_d0.isValid()) return;
    
    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByToken(tok_generalTrk_, tracks);
    
    for(int i=0;i<ptbin_n_;i++)
    {
        pVect_trg_d0[i] = new vector<TVector3>;
        pVect_dau_d0[i] = new vector<TVector3>;
        pVect_trg_d0_bkg[i] = new vector<TVector3>;
        pVect_dau_d0_bkg[i] = new vector<TVector3>;
     }
    
    pVect_ass = new vector<TVector3>;
    

    //track selection
    int nMult_ass_good = 0;
    for(unsigned it=0; it<tracks->size(); ++it){
        
        const reco::Track & trk = (*tracks)[it];
        
        math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
        
        double dzvtx = trk.dz(bestvtx);
        double dxyvtx = trk.dxy(bestvtx);
        double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
        double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError);
        
        if(!trk.quality(reco::TrackBase::highPurity)) continue;
        if(fabs(trk.ptError())/trk.pt()>0.10) continue;
        if(fabs(dzvtx/dzerror) > 3) continue;
        if(fabs(dxyvtx/dxyerror) > 3) continue;
        
        double eta = trk.eta();
        double pt  = trk.pt();
        
        if(fabs(eta)>2.4) continue;
        if(pt<=0.4) continue;
        nMult_ass_good++;
    }
    hMult->Fill(nMult_ass_good);

    if(nMult_ass_good<multMax_ && nMult_ass_good>=multMin_){
        for(unsigned it=0; it<v0candidates_d0->size(); ++it){

            const reco::VertexCompositeCandidate & trk = (*v0candidates_d0)[it];
            
            double secvz=-999.9, secvx=-999.9, secvy=-999.9;
            //double secvzError=-999.9, secvxError=-999.9, secvyError=-999.9;
            
            secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();
            //secvzError = trk.vertexCovariance(2,2); secvxError = trk.vertexCovariance(0,0); secvyError = trk.vertexCovariance(1,1);
            
            double eta = trk.eta();
            double phi = trk.phi();
            double pt  = trk.pt();
            double mass = trk.mass();
            double px = trk.px();
            double py = trk.py();
            double pz = trk.pz();
            double rapidity = trk.rapidity();
            
            if(fabs(rapidity)>d0_rapidity_) continue;
            
            TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
            TVector3 secvec(px,py,pz);
            
            double agl = secvec.Angle(ptosvec);
            
            //if(agl<=0.999) continue;
            
            typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
            typedef ROOT::Math::SVector<double, 3> SVector3;
            
            SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
            SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);
            
            double dl = ROOT::Math::Mag(distanceVector);
            double dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl;
            
            double dlos = dl/dlerror;
            
            //if(dlos<=5) continue;
            
            double VtxProb = TMath::Prob(trk.vertexChi2(),trk.vertexNdof());
            
            const reco::Candidate * dau1 = trk.daughter(0);
            const reco::Candidate * dau2 = trk.daughter(1);
            
            auto d1 = dau1->get<reco::TrackRef>();
            auto d2 = dau2->get<reco::TrackRef>();
            
            double pxd1 = dau1->px();
            double pyd1 = dau1->py();
            double pzd1 = dau1->pz();
            double pd1 = dau1->p();
            double pxd2 = dau2->px();
            double pyd2 = dau2->py();
            double pzd2 = dau2->pz();
            double pd2 = dau2->p();
            
            TVector3 dauvec1(pxd1,pyd1,pzd1);
            TVector3 dauvec2(pxd2,pyd2,pzd2);
            
            //efficiency
            //double effks = effhisto_d0->GetBinContent(effhisto_d0->FindBin(eta,pt));
            double effks = 1.0;

            double eta_dau1 = dau1->eta();
            double phi_dau1 = dau1->phi();
            double pt_dau1 = dau1->pt();
            double pt_err_dau1 = d1->ptError();
            int nhit_dau1 = d1->numberOfValidHits();
            
            double eta_dau2 = dau2->eta();
            double phi_dau2 = dau2->phi();
            double pt_dau2 = dau2->pt();
            double pt_err_dau2 = d2->ptError();
            int nhit_dau2 = d2->numberOfValidHits();
            
            TVector3 pvector;
            pvector.SetPtEtaPhi(pt,eta,phi);
            
            TVector3 pvector_dau1;
            pvector_dau1.SetPtEtaPhi(pt_dau1,eta_dau1,phi_dau1);
            
            TVector3 pvector_dau2;
            pvector_dau2.SetPtEtaPhi(pt_dau2,eta_dau2,phi_dau2);
            
            for(int i=0;i<ptbin_n_;i++)
            {
                if(fabs(eta)<=d0_eta_ && pt<=ptcut_d0_[i+1] && pt>=ptcut_d0_[i] && agl<angle_d0_[i] && dlos>dls_d0_[i] && VtxProb>vtxprob_d0_[i] && nhit_dau1>=d0_dau_nhit_ && nhit_dau2>=d0_dau_nhit_ && fabs(eta_dau1)<d0_dau_eta_ && fabs(eta_dau2)<d0_dau_eta_ && pt_dau1>daupt_d0_[i] && pt_dau2>daupt_d0_[i] && fabs(pt_err_dau1)<0.1 && fabs(pt_err_dau2)<0.1){//FIXME pTerr/pT
                    hMass_d0[i]->Fill(mass);
                    if(mass<=(mean_d0_[i]+peakFactor_*sigma_d0_[i]) && mass>=(mean_d0_[i]-peakFactor_*sigma_d0_[i])){
                        pVect_trg_d0[i]->push_back(pvector);
                        pVect_dau_d0[i]->push_back(pvector_dau1);
                        pVect_dau_d0[i]->push_back(pvector_dau2);
                        hPt_d0[i]->Fill(pt,1.0/effks);
                        hEta_d0[i]->Fill(eta,1.0/effks);
                        double KET = sqrt(mass*mass + pt*pt) - mass;
                        hKET_d0[i]->Fill(KET,1.0/effks);
                    }
                    if(mass<=(mean_d0_[i]-sideFactor_*sigma_d0_[i]) || mass>=(mean_d0_[i]+sideFactor_*sigma_d0_[i])){
                        pVect_trg_d0_bkg[i]->push_back(pvector);
                        pVect_dau_d0_bkg[i]->push_back(pvector_dau1);
                        pVect_dau_d0_bkg[i]->push_back(pvector_dau2);
                        hPt_d0_bkg[i]->Fill(pt,1.0/effks);
                        hEta_d0_bkg[i]->Fill(eta,1.0/effks);
                        double KET = sqrt(mass*mass + pt*pt) - mass;
                        hKET_d0_bkg[i]->Fill(KET,1.0/effks);
                    }
                }
            }
        }

        for(unsigned it=0; it<tracks->size(); ++it){
            
            const reco::Track & trk = (*tracks)[it];
            
            math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
            
            double dzvtx = trk.dz(bestvtx);
            double dxyvtx = trk.dxy(bestvtx);
            double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
            double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError);
            
            if(!trk.quality(reco::TrackBase::highPurity)) continue;
            if(fabs(trk.ptError())/trk.pt()>0.10) continue;
            if(fabs(dzvtx/dzerror) > 3) continue;
            if(fabs(dxyvtx/dxyerror) > 3) continue;
            
            
            double eta = trk.eta();
            double phi = trk.phi();
            double pt  = trk.pt();
            
            TVector3 pvector;
            pvector.SetPtEtaPhi(pt,eta,phi);
            if(eta<=etaMax_ass_ && eta>=etaMin_ass_ && pt<=ptMax_ass_ && pt>=ptMin_ass_) pVect_ass->push_back(pvector);
        }

        //Calculating signal
        int nMult_ass = (int)pVect_ass->size();
        hMult_ass->Fill(nMult_ass);

        for(int i=0; i<ptbin_n_; i++)
        {
            int nMult_trg_d0 = (int)pVect_trg_d0[i]->size();
            double nMult_trg_eff_d0=0;
            
            for(int ntrg=0;ntrg<nMult_trg_d0;ntrg++)
            {
                TVector3 pvector_trg = (*pVect_trg_d0[i])[ntrg];
                double pt_trg = pvector_trg.Pt();
                double eta_trg = pvector_trg.Eta();
                
                //double effks = effhisto_d0->GetBinContent(effhisto_d0->FindBin(eta_trg,pt_trg));
                double effks = 1.0;

                nMult_trg_eff_d0 = nMult_trg_eff_d0 + 1.0/effks;
            }
            
            hMult_d0[i]->Fill(nMult_trg_d0);
            
            for(int ntrg=0;ntrg<nMult_trg_d0;ntrg++)
            {
                TVector3 pvector_trg = (*pVect_trg_d0[i])[ntrg];
                double eta_trg = pvector_trg.Eta();
                double phi_trg = pvector_trg.Phi();
                double pt_trg = pvector_trg.Pt();
                
                //double effks = effhisto_d0->GetBinContent(effhisto_d0->FindBin(eta_trg,pt_trg));
                double effks = 1.0;

                TVector3 pvector_trg_dau1 = (*pVect_dau_d0[i])[2*ntrg];
                double eta_trg_dau1 = pvector_trg_dau1.Eta();
                double phi_trg_dau1 = pvector_trg_dau1.Phi();
                
                TVector3 pvector_trg_dau2 = (*pVect_dau_d0[i])[2*ntrg+1];
                double eta_trg_dau2 = pvector_trg_dau2.Eta();
                double phi_trg_dau2 = pvector_trg_dau2.Phi();
                
                for(int nass=0;nass<nMult_ass;nass++)
                {
                    TVector3 pvector_ass = (*pVect_ass)[nass];
                    double eta_ass = pvector_ass.Eta();
                    double phi_ass = pvector_ass.Phi();
                    double pt_ass = pvector_ass.Pt();

                    double effweight_ass = effhisto->GetBinContent(effhisto->FindBin(eta_ass,pt_ass));
                    
                    if(rejectDaughter_)
                    {
                    if(fabs(eta_ass-eta_trg_dau1)<0.03 && fabs(phi_ass-phi_trg_dau1)<0.03) continue;
                    if(fabs(eta_ass-eta_trg_dau2)<0.03 && fabs(phi_ass-phi_trg_dau2)<0.03) continue;
                    }
   
                    double deltaEta=eta_ass-eta_trg;
                    double deltaPhi=phi_ass-phi_trg;
                    if(deltaPhi>PI)
                        deltaPhi=deltaPhi-2*PI;
                    if(deltaPhi<-PI)
                        deltaPhi=deltaPhi+2*PI;
                    if(deltaPhi>-PI && deltaPhi<-PI/2.)
                        deltaPhi=deltaPhi+2*PI;
                    
                    //if(deltaEta==0 && deltaPhi==0) continue;
                    hSignal_d0[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_eff_d0/effks/effweight_ass);
                }
            }
        }

        for(int i=0; i<ptbin_n_; i++)
        {
            int nMult_trg_d0 = (int)pVect_trg_d0_bkg[i]->size();
            
            double nMult_trg_eff_d0=0;
            
            for(int ntrg=0;ntrg<nMult_trg_d0;ntrg++)
            {
                TVector3 pvector_trg = (*pVect_trg_d0_bkg[i])[ntrg];
                double pt_trg = pvector_trg.Pt();
                double eta_trg = pvector_trg.Eta();
                
                //double effks = effhisto_d0->GetBinContent(effhisto_d0->FindBin(eta_trg,pt_trg));
                double effks = 1.0;

                nMult_trg_eff_d0 = nMult_trg_eff_d0 + 1.0/effks;
            }
            
            hMult_d0_bkg[i]->Fill(nMult_trg_d0);
            
            for(int ntrg=0;ntrg<nMult_trg_d0;ntrg++)
            {
                TVector3 pvector_trg = (*pVect_trg_d0_bkg[i])[ntrg];
                double eta_trg = pvector_trg.Eta();
                double phi_trg = pvector_trg.Phi();
                double pt_trg = pvector_trg.Pt();

                //double effks = effhisto_d0->GetBinContent(effhisto_d0->FindBin(eta_trg,pt_trg));
                double effks = 1.0;

                TVector3 pvector_trg_dau1 = (*pVect_dau_d0_bkg[i])[2*ntrg];
                double eta_trg_dau1 = pvector_trg_dau1.Eta();
                double phi_trg_dau1 = pvector_trg_dau1.Phi();
                
                TVector3 pvector_trg_dau2 = (*pVect_dau_d0_bkg[i])[2*ntrg+1];
                double eta_trg_dau2 = pvector_trg_dau2.Eta();
                double phi_trg_dau2 = pvector_trg_dau2.Phi();

                for(int nass=0;nass<nMult_ass;nass++)
                {
                    TVector3 pvector_ass = (*pVect_ass)[nass];
                    double eta_ass = pvector_ass.Eta();
                    double phi_ass = pvector_ass.Phi();
                    double pt_ass = pvector_ass.Pt();
                    double effweight_ass = effhisto->GetBinContent(effhisto->FindBin(eta_ass,pt_ass));

                    if(rejectDaughter_)
                    {
                    if(fabs(eta_ass-eta_trg_dau1)<0.03 && fabs(phi_ass-phi_trg_dau1)<0.03) continue;
                    if(fabs(eta_ass-eta_trg_dau2)<0.03 && fabs(phi_ass-phi_trg_dau2)<0.03) continue;
                    }
   
                    double deltaEta=eta_ass-eta_trg;
                    double deltaPhi=phi_ass-phi_trg;
                    if(deltaPhi>PI)
                        deltaPhi=deltaPhi-2*PI;
                    if(deltaPhi<-PI)
                        deltaPhi=deltaPhi+2*PI;
                    if(deltaPhi>-PI && deltaPhi<-PI/2.)
                        deltaPhi=deltaPhi+2*PI;

                    //if(deltaEta==0 && deltaPhi==0) continue;
                    hSignal_d0_bkg[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_eff_d0/effks/effweight_ass);
                }
            }
        }
        
        for(int i=0; i<ptbin_n_; i++)
        {
            pVectVect_trg_d0[i]->push_back(*pVect_trg_d0[i]);
            pVectVect_dau_d0[i]->push_back(*pVect_dau_d0[i]);
            delete pVect_trg_d0[i];
            delete pVect_dau_d0[i];
            pVectVect_trg_d0_bkg[i]->push_back(*pVect_trg_d0_bkg[i]);
            pVectVect_dau_d0_bkg[i]->push_back(*pVect_dau_d0_bkg[i]);
            delete pVect_trg_d0_bkg[i];
            delete pVect_dau_d0_bkg[i];
        }
        
        pVectVect_ass->push_back(*pVect_ass);
        zvtxVect->push_back(bestvz);
        delete pVect_ass;
    }
}


// ------------ method called once each job just before starting event
//loop  ------------
void
D0Corr::beginJob()
{
    edm::Service<TFileService> fs;
    
    TH1::SetDefaultSumw2();

    edm::FileInPath fip("D0Corr/D0Corr/data/Hijing_8TeV_dataBS.root");
    TFile f(fip.fullPath().c_str(),"READ");
    effhisto = (TH2F*)f.Get("rTotalEff3D_0");
    
    edm::FileInPath fip1("D0Corr/D0Corr/data/Hijing_8TeV_dataBS.root");
    TFile f1(fip1.fullPath().c_str(),"READ");
    effhisto_d0 = (TH2D*)f1.Get("rTotalEff3D_0");
    
    hMult = fs->make<TH1D>("mult",";N",600,0,600);
    hMult_ass = fs->make<TH1D>("mult_ass",";N",600,0,600);
    
    for(int i=0; i<ptbin_n_; i++)
    {
        hKET_d0[i] = fs->make<TH1D>(Form("KETD0_pt%d",i),";GeV",2500,0,12.5);
        hPt_d0[i] = fs->make<TH1D>(Form("PtD0_pt%d",i),";GeV",2500,0,12.5);
        hEta_d0[i] = fs->make<TH1D>(Form("EtaD0_pt%d",i),";eta",24,-2.4,2.4);
        hMass_d0[i] = fs->make<TH1D>(Form("massD0_pt%d",i),";GeV",200,1.5,2.5);
        hMult_d0[i] = fs->make<TH1D>(Form("mult_d0_pt%d",i),";N",250,0,250);
        hSignal_d0[i] = fs->make<TH2D>(Form("signalD0_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        hBackground_d0[i] = fs->make<TH2D>(Form("backgroundD0_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        pVectVect_trg_d0[i] = new vector< vector<TVector3> >;
        pVectVect_dau_d0[i] = new vector< vector<TVector3> >;
        
        //bkg left
        hKET_d0_bkg[i] = fs->make<TH1D>(Form("KETD0_bkg_pt%d",i),";GeV",2500,0,12.5);
        hPt_d0_bkg[i] = fs->make<TH1D>(Form("PtD0_bkg_pt%d",i),";GeV",2500,0,12.5);
        hEta_d0_bkg[i] = fs->make<TH1D>(Form("EtaD0_bkg_pt%d",i),";GeV",24,-2.4,2.4);
        hMult_d0_bkg[i] = fs->make<TH1D>(Form("mult_d0_bkg_pt%d",i),";N",250,0,250);
        hSignal_d0_bkg[i] = fs->make<TH2D>(Form("signalD0_bkg__pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        hBackground_d0_bkg[i] = fs->make<TH2D>(Form("backgroundD0_bkg_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        pVectVect_trg_d0_bkg[i] = new vector< vector<TVector3> >;
        pVectVect_dau_d0_bkg[i] = new vector< vector<TVector3> >;
    }
    
    pVectVect_ass = new vector< vector<TVector3> >;
    zvtxVect = new vector<double>;
    
}

// ------------ method called once each job just after ending the event
//loop  ------------
void
D0Corr::endJob() {
    //Calculating background
    int nevttotal_ass = (int)pVectVect_ass->size();
    
    for(int i=0;i<ptbin_n_;i++)
    {
        int nevttotal_trg_d0 = (int)pVectVect_trg_d0[i]->size();
        
        for(int nround=0;nround<bkgFactor_;nround++)
        {
            int ncount = 0;
            for(int nevt_trg=0; nevt_trg<nevttotal_trg_d0; nevt_trg++)
            {
                int nevt_ass = gRandom->Integer(nevttotal_ass);
                if(nevt_trg == nevt_ass) { nevt_trg--; continue; }
                if(fabs((*zvtxVect)[nevt_trg]-(*zvtxVect)[nevt_ass])>0.5) {
                    nevt_trg--;
                    ncount++;
                    if(ncount>5000) {nevt_trg++; ncount = 0;}
                    continue; }
                
                vector<TVector3> pVectTmp_trg = (*pVectVect_trg_d0[i])[nevt_trg];
                vector<TVector3> pVectTmp_dau = (*pVectVect_dau_d0[i])[nevt_trg];
                vector<TVector3> pVectTmp_ass = (*pVectVect_ass)[nevt_ass];
                int nMult_trg = pVectTmp_trg.size();
                int nMult_ass = pVectTmp_ass.size();
                
                double nMult_trg_eff=0;
                
                for(int ntrg=0;ntrg<nMult_trg;ntrg++)
                {
                    TVector3 pvectorTmp_trg = pVectTmp_trg[ntrg];
                    double pt_trg = pvectorTmp_trg.Pt();
                    double eta_trg = pvectorTmp_trg.Eta();
                    
                    //double effks = effhisto_d0->GetBinContent(effhisto_d0->FindBin(eta_trg,pt_trg));
                    double effks = 1.0;

                    nMult_trg_eff = nMult_trg_eff + 1.0/effks;
                }
                
                for(int ntrg=0;ntrg<nMult_trg;ntrg++)
                {
                    TVector3 pvectorTmp_trg = pVectTmp_trg[ntrg];
                    double eta_trg = pvectorTmp_trg.Eta();
                    double phi_trg = pvectorTmp_trg.Phi();
                    double pt_trg = pvectorTmp_trg.Pt();
                    
                    //double effks = effhisto_d0->GetBinContent(effhisto_d0->FindBin(eta_trg,pt_trg));
                    double effks = 1.0;

                    TVector3 pvector_trg_dau1 = pVectTmp_dau[2*ntrg];
                    double eta_trg_dau1 = pvector_trg_dau1.Eta();
                    double phi_trg_dau1 = pvector_trg_dau1.Phi();
                    
                    TVector3 pvector_trg_dau2 = pVectTmp_dau[2*ntrg+1];
                    double eta_trg_dau2 = pvector_trg_dau2.Eta();
                    double phi_trg_dau2 = pvector_trg_dau2.Phi();

                    for(int nass=0;nass<nMult_ass;nass++)
                    {
                        TVector3 pvectorTmp_ass = pVectTmp_ass[nass];
                        double eta_ass = pvectorTmp_ass.Eta();
                        double phi_ass = pvectorTmp_ass.Phi();
                        double pt_ass = pvectorTmp_ass.Pt();
                        
                        double effweight_ass = effhisto->GetBinContent(effhisto->FindBin(eta_ass,pt_ass));
                        
                        if(rejectDaughter_)
                        {
                        if(fabs(eta_ass - eta_trg_dau1)<0.03 && fabs(phi_ass - phi_trg_dau1)<0.03) continue;
                        if(fabs(eta_ass - eta_trg_dau2)<0.03 && fabs(phi_ass - phi_trg_dau2)<0.03) continue;
                        }

                        double deltaEta=eta_ass-eta_trg;
                        double deltaPhi=phi_ass-phi_trg;
                        if(deltaPhi>PI)
                            deltaPhi=deltaPhi-2*PI;
                        if(deltaPhi<-PI)
                            deltaPhi=deltaPhi+2*PI;
                        if(deltaPhi>-PI && deltaPhi<-PI/2.)
                            deltaPhi=deltaPhi+2*PI;
                        
                        //if(deltaEta==0 && deltaPhi==0) continue;
                        hBackground_d0[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_eff/effks/effweight_ass);
                    }
                }
            }
        }
    }
    
    for(int i=0;i<ptbin_n_;i++)
    {
        int nevttotal_trg_d0 = (int)pVectVect_trg_d0_bkg[i]->size();
        
        for(int nround=0;nround<bkgFactor_;nround++)
        {
            int ncount = 0;
            for(int nevt_trg=0; nevt_trg<nevttotal_trg_d0; nevt_trg++)
            {
                int nevt_ass = gRandom->Integer(nevttotal_ass);
                if(nevt_trg == nevt_ass) { nevt_trg--; continue; }
                if(fabs((*zvtxVect)[nevt_trg]-(*zvtxVect)[nevt_ass])>0.5) {
                    nevt_trg--;
                    ncount++;
                    if(ncount>5000) {nevt_trg++; ncount = 0;}
                    continue; }
                
                vector<TVector3> pVectTmp_trg = (*pVectVect_trg_d0_bkg[i])[nevt_trg];
                vector<TVector3> pVectTmp_dau = (*pVectVect_dau_d0_bkg[i])[nevt_trg];
                vector<TVector3> pVectTmp_ass = (*pVectVect_ass)[nevt_ass];
                int nMult_trg = pVectTmp_trg.size();
                int nMult_ass = pVectTmp_ass.size();
                
                double nMult_trg_eff=0;
                
                for(int ntrg=0;ntrg<nMult_trg;ntrg++)
                {
                    TVector3 pvectorTmp_trg = pVectTmp_trg[ntrg];
                    double pt_trg = pvectorTmp_trg.Pt();
                    double eta_trg = pvectorTmp_trg.Eta();
                    
                    //double effks = effhisto_d0->GetBinContent(effhisto_d0->FindBin(eta_trg,pt_trg));
                    double effks = 1.0;
                    nMult_trg_eff = nMult_trg_eff + 1.0/effks;
                }
                
                for(int ntrg=0;ntrg<nMult_trg;ntrg++)
                {
                    TVector3 pvectorTmp_trg = pVectTmp_trg[ntrg];
                    double eta_trg = pvectorTmp_trg.Eta();
                    double phi_trg = pvectorTmp_trg.Phi();
                    double pt_trg = pvectorTmp_trg.Pt();
                    
                    //double effks = effhisto_d0->GetBinContent(effhisto_d0->FindBin(eta_trg,pt_trg));
                    double effks = 1.0;

                    TVector3 pvector_trg_dau1 = pVectTmp_dau[2*ntrg];
                    double eta_trg_dau1 = pvector_trg_dau1.Eta();
                    double phi_trg_dau1 = pvector_trg_dau1.Phi();
                    
                    TVector3 pvector_trg_dau2 = pVectTmp_dau[2*ntrg+1];
                    double eta_trg_dau2 = pvector_trg_dau2.Eta();
                    double phi_trg_dau2 = pvector_trg_dau2.Phi();
                    
                    for(int nass=0;nass<nMult_ass;nass++)
                    {
                        TVector3 pvectorTmp_ass = pVectTmp_ass[nass];
                        double eta_ass = pvectorTmp_ass.Eta();
                        double phi_ass = pvectorTmp_ass.Phi();
                        double pt_ass = pvectorTmp_ass.Pt();
                        
                        double effweight_ass = effhisto->GetBinContent(effhisto->FindBin(eta_ass,pt_ass));

                        if(rejectDaughter_)
                        {
                        if(fabs(eta_ass - eta_trg_dau1)<0.03 && fabs(phi_ass - phi_trg_dau1)<0.03) continue;
                        if(fabs(eta_ass - eta_trg_dau2)<0.03 && fabs(phi_ass - phi_trg_dau2)<0.03) continue;
                        }

                        double deltaEta=eta_ass-eta_trg;
                        double deltaPhi=phi_ass-phi_trg;
                        if(deltaPhi>PI)
                            deltaPhi=deltaPhi-2*PI;
                        if(deltaPhi<-PI)
                            deltaPhi=deltaPhi+2*PI;
                        if(deltaPhi>-PI && deltaPhi<-PI/2.)
                            deltaPhi=deltaPhi+2*PI;
                        
                        //if(deltaEta==0 && deltaPhi==0) continue;
                        hBackground_d0_bkg[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_eff/effks/effweight_ass);
                    }
                }
            }
        }
    }
}

//define this as a plug-in
DEFINE_FWK_MODULE(D0Corr);



