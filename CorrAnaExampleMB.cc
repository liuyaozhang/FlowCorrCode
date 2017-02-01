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

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
//
// class declaration
//

#define PI 3.1416

using namespace std;

class CorrAnaExampleMB : public edm::EDAnalyzer {
public:
    explicit CorrAnaExampleMB(const edm::ParameterSet&);
    ~CorrAnaExampleMB();
    
    
private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;
    
    // ----------member data ---------------------------
    TH1D* hMult_selected;
    TH1D* hMult;
    TH1D* hMult_good;
    TH1D* hMult_assoc;
    TH1D* hPt;
    TH1D* hEta;
    TH2D* hSignal;
    TH2D* hBackground;
    //TH2D* hCorrelation;
    
    //TH1D* hMass_ks[9];
    //TH1D* hMass_la[9];
    
    TH2F* effhisto;
    TH1D* hNvtx;
    double ptcut[10];
    bool isMC_;
    
    vector<TVector3> *pVect_trg;
    vector< vector<TVector3> > *pVectVect_trg;
    vector<TVector3> *pVect_ass;
    vector< vector<TVector3> > *pVectVect_ass;
    vector<double> *zvtxVect;
    
    double etaMin_trg_;
    double etaMax_trg_;
    double etaMin_ass_;
    double etaMax_ass_;
    double ptMin_trg_;
    double ptMax_trg_;
    double ptMin_ass_;
    double ptMax_ass_;
    int bkgFactor_;
    double multMax_;
    double multMin_;
    double PUfrac_;
    double dzCut_;
    double toppt_min_;
    double toppt_max_;
    edm::Service<TFileService> fs;
    
    edm::EDGetTokenT<reco::VertexCollection> tok_offlinePV_;
    edm::EDGetTokenT<reco::TrackCollection> tok_generalTrk_;
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
CorrAnaExampleMB::CorrAnaExampleMB(const edm::ParameterSet& iConfig)

{
    //now do what ever initialization is needed
    etaMin_trg_ = iConfig.getUntrackedParameter<double>("etaMin_trg", -2.4);
    etaMax_trg_ = iConfig.getUntrackedParameter<double>("etaMax_trg", 2.4);
    etaMin_ass_ = iConfig.getUntrackedParameter<double>("etaMin_ass", -2.4);
    etaMax_ass_ = iConfig.getUntrackedParameter<double>("etaMax_ass", 2.4);
    ptMin_trg_ = iConfig.getUntrackedParameter<double>("ptMin_trg", 0.3);
    ptMax_trg_ = iConfig.getUntrackedParameter<double>("ptMax_trg", 3.0);
    ptMin_ass_ = iConfig.getUntrackedParameter<double>("ptMin_ass", 0.3);
    ptMax_ass_ = iConfig.getUntrackedParameter<double>("ptMax_ass", 3.0);
    bkgFactor_ = iConfig.getUntrackedParameter<int>("bkgFactor", 10);
    multMax_ = iConfig.getUntrackedParameter<double>("multMax", 220);
    multMin_ = iConfig.getUntrackedParameter<double>("multMin", 185);
    PUfrac_ = iConfig.getUntrackedParameter<double>("PUfrac",0.1);
    dzCut_ = iConfig.getUntrackedParameter<double>("dzCut",9999.0);
    isMC_ = iConfig.getUntrackedParameter<bool>("isMC");
    toppt_min_ = iConfig.getUntrackedParameter<double>("toppt_min",0.0);
    toppt_max_ = iConfig.getUntrackedParameter<double>("toppt_max",9999.0);
    
    tok_offlinePV_ = consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"));
    tok_generalTrk_ = consumes<reco::TrackCollection>(edm::InputTag("generalTracks"));
}


CorrAnaExampleMB::~CorrAnaExampleMB()
{
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
CorrAnaExampleMB::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    
    //double ptcut[10] = {0.2,0.6,1.2,1.8,2.6,3.6,4.6,6.0,9.0,12.0};
    
    if(isMC_)
    {
        edm::Handle<edm::HepMCProduct> hepmc;
        iEvent.getByLabel("generator",hepmc);
        
        const HepMC::GenEvent *myGenEvent = hepmc->GetEvent();
        int procID = myGenEvent->signal_process_id();
        
        if(procID==92 || procID==93 || procID==94) return;
    }
    
    // select on requirement of valid vertex
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(tok_offlinePV_,vertices);
    double bestvz=-999.9, bestvx=-999.9, bestvy=-999.9;
    double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
    
    //reject PU
    hNvtx->Fill(vertices->size());
    
    int nTrack1 = 0;
    int nTrack2 = 0;
    
    for(unsigned it=0; it<vertices->size(); ++it)
    {
        const reco::Vertex & vtx1 = (*vertices)[it];
        int nTrack = vtx1.tracksSize();
        if(nTrack>nTrack1)
        {
            nTrack2 = nTrack1;
            nTrack1 = nTrack;
        }
        if(nTrack<nTrack1 && nTrack>nTrack2)
        {
            nTrack2 = nTrack;
        }
    }
    
    //if(!(vertices->size()==1 || nTrack1*PUfrac_ > nTrack2)) return;
    
    
    const reco::Vertex & vtx = (*vertices)[0];
    bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
    bestvzError = vtx.zError(); bestvxError = vtx.xError(); bestvyError = vtx.yError();
    
    if(bestvz < -15.0 || bestvz>15.0) return;
    //if(bestvz < -10.4 || bestvz>9.6) return;
    
    /*edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates_ks;
     iEvent.getByLabel("generalV0CandidatesNew","Kshort",v0candidates_ks);
     if(!v0candidates_ks.isValid()) return;
     
     edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates_la;
     iEvent.getByLabel("generalV0CandidatesNew","Lambda",v0candidates_la);
     if(!v0candidates_la.isValid()) return;
     */
    pVect_trg = new vector<TVector3>;
    pVect_ass = new vector<TVector3>;
    //----- loop over tracks -----
    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByToken(tok_generalTrk_, tracks);
    
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
        //if(fabs(dzvtx) > dzCut_) continue;
        
        double eta = trk.eta();
        double pt  = trk.pt();
        
        if(fabs(eta)>2.4) continue;
        if(pt<=0.4) continue;
        nMult_ass_good++;
    }
    hMult_selected->Fill(nMult_ass_good);
    
    double toppt=-999;
    
    if(nMult_ass_good<multMax_ && nMult_ass_good>=multMin_){
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
            
            if(eta<=etaMax_trg_ && eta>=etaMin_trg_ && pt>toppt) toppt = pt;
        }
    }
    
    if(toppt < toppt_min_ || toppt > toppt_max_) return;
    
    if(nMult_ass_good<multMax_ && nMult_ass_good>=multMin_){
        hMult_good->Fill(nMult_ass_good);
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
            if(fabs(dzvtx) > dzCut_) continue;
            
            double eta = trk.eta();
            double phi = trk.phi();
            double pt  = trk.pt();
            
            TVector3 pvector;
            pvector.SetPtEtaPhi(pt,eta,phi);
            
            double effweight = effhisto->GetBinContent(effhisto->FindBin(eta,pt));
            
            if(eta<=etaMax_trg_ && eta>=etaMin_trg_ && pt<=ptMax_trg_ && pt>=ptMin_trg_){ pVect_trg->push_back(pvector); hPt->Fill(pt,1.0/effweight); hEta->Fill(eta,1.0/effweight);}
            if(eta<=etaMax_ass_ && eta>=etaMin_ass_ && pt<=ptMax_ass_ && pt>=ptMin_ass_) pVect_ass->push_back(pvector);
        }
        // Calculating signal
        int nMult_trg = (int)pVect_trg->size();
        int nMult_ass = (int)pVect_ass->size();
        hMult->Fill(nMult_trg);
        
        double nMult_trg_eff=0;
        
        for(int ntrg=0;ntrg<nMult_trg;ntrg++)
        {
            TVector3 pvector_trg = (*pVect_trg)[ntrg];
            double eta_trg = pvector_trg.Eta();
            //double phi_trg = pvector_trg.Phi();
            double pt_trg = pvector_trg.Pt();
            
            double effweight_trg = effhisto->GetBinContent(effhisto->FindBin(eta_trg,pt_trg));
            
            nMult_trg_eff = nMult_trg_eff + 1.0/effweight_trg;
        }
        
        
        for(int ntrg=0;ntrg<nMult_trg;ntrg++)
        {
            int nMult_ass_pair = (int)pVect_ass->size();
            TVector3 pvector_trg = (*pVect_trg)[ntrg];
            double eta_trg = pvector_trg.Eta();
            double phi_trg = pvector_trg.Phi();
            double pt_trg = pvector_trg.Pt();
            
            double effweight_trg = effhisto->GetBinContent(effhisto->FindBin(eta_trg,pt_trg));
            
            for(int nass=0;nass<nMult_ass;nass++)
            {
                TVector3 pvector_ass = (*pVect_ass)[nass];
                double eta_ass = pvector_ass.Eta();
                double phi_ass = pvector_ass.Phi();
                double pt_ass = pvector_ass.Pt();
                
                double effweight_ass = effhisto->GetBinContent(effhisto->FindBin(eta_ass,pt_ass));
                
                double deltaEta=eta_ass-eta_trg;
                double deltaPhi=phi_ass-phi_trg;
                if(deltaPhi>PI)
                    deltaPhi=deltaPhi-2*PI;
                if(deltaPhi<-PI)
                    deltaPhi=deltaPhi+2*PI;
                if(deltaPhi>-PI && deltaPhi<-PI/2.)
                    deltaPhi=deltaPhi+2*PI;
                
                //if(deltaEta==0 && deltaPhi==0) continue;
                if(fabs(deltaEta)==0 && fabs(deltaPhi)==0){
                    nMult_ass_pair = nMult_ass_pair - 1;
                    continue;
                }
                hSignal->Fill(deltaEta,deltaPhi,1.0/nMult_trg_eff/effweight_trg/effweight_ass);
            }
            hMult_assoc->Fill(nMult_ass_pair);
        }
        pVectVect_trg->push_back(*pVect_trg);
        pVectVect_ass->push_back(*pVect_ass);
        zvtxVect->push_back(bestvz);
        
        delete pVect_trg;
        delete pVect_ass;
    }
}

// ------------ method called once each job just before starting event loop  ------------
void
CorrAnaExampleMB::beginJob(){
    
    TH1D::SetDefaultSumw2();
    
    edm::FileInPath fip("pPbCorr/pPbCorr/data/Hydjet_ppReco_dataBS.root");
    TFile f(fip.fullPath().c_str(),"READ");
    effhisto = (TH2F*)f.Get("rTotalEff3D_0");
    
    hMult_selected = fs->make<TH1D>("mult_selected",";N",300,0,300);
    hMult_good = fs->make<TH1D>("mult_good",";N",300,0,300);
    hMult_assoc = fs->make<TH1D>("mult_assoc",";N",300,0,300);
    hMult = fs->make<TH1D>("mult",";N",1500,0,1500);
    hSignal = fs->make<TH2D>("signal",";#Delta#eta;#Delta#phi",
                             33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
    hBackground = fs->make<TH2D>("background",";#Delta#eta;#Delta#phi",
                                 33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
    hNvtx = fs->make<TH1D>("mult_vtx",";N",20,0,20);
    hPt = fs->make<TH1D>("pT",";pT",20000,0,10);
    hEta = fs->make<TH1D>("eta",";eta",24,-2.4,2.4);
    
    //hCorrelation = fs->make<TH2D>("correlation",";#Delta#eta;#Delta#phi",
    //                              27,-4.05,4.05,33,-(0.5+1/32)*PI,(1.5+1/32)*PI);
    /*for(int i=0; i<9; i++)
     {
     hMass_ks[i] = fs->make<TH1D>(Form("masskshort_pt%d",i),";GeV",2000,0,1.0);
     hMass_la[i] = fs->make<TH1D>(Form("masslambda_pt%d",i),";GeV",2000,0.5,1.5);
     }*/
    pVectVect_trg = new vector< vector<TVector3> >;
    pVectVect_ass = new vector< vector<TVector3> >;
    zvtxVect = new vector<double>;
}

// ------------ method called once each job just after ending the event loop  ------------
void
CorrAnaExampleMB::endJob() {
    // Calculate background
    int nevttotal_trg = (int)pVectVect_trg->size();
    int nevttotal_ass = (int)pVectVect_ass->size();
    
    for(int nround=0;nround<bkgFactor_;nround++)
    {
        int ncount = 0;
        for(int nevt_ass=0; nevt_ass<nevttotal_ass; nevt_ass++)
        {
            int nevt_trg = gRandom->Integer(nevttotal_trg);
            if(nevt_trg == nevt_ass) { nevt_ass--; continue; }
            if(fabs((*zvtxVect)[nevt_trg]-(*zvtxVect)[nevt_ass])>0.5) {
                nevt_ass--;
                ncount++;
                if(ncount>5000) {nevt_ass++; ncount = 0;}
                continue; }
            
            vector<TVector3> pVectTmp_trg = (*pVectVect_trg)[nevt_trg];
            vector<TVector3> pVectTmp_ass = (*pVectVect_ass)[nevt_ass];
            int nMult_trg = pVectTmp_trg.size();
            int nMult_ass = pVectTmp_ass.size();
            
            double nMult_trg_eff=0;
            
            for(int ntrg=0;ntrg<nMult_trg;ntrg++)
            {
                TVector3 pvectorTmp_trg = pVectTmp_trg[ntrg];
                double eta_trg = pvectorTmp_trg.Eta();
                //double phi_trg = pvectorTmp_trg.Phi();
                double pt_trg = pvectorTmp_trg.Pt();
                
                double effweight_trg = effhisto->GetBinContent(effhisto->FindBin(eta_trg,pt_trg));
                
                nMult_trg_eff = nMult_trg_eff + 1.0/effweight_trg;
            }
            
            for(int ntrg=0;ntrg<nMult_trg;ntrg++)
            {
                TVector3 pvectorTmp_trg = pVectTmp_trg[ntrg];
                double eta_trg = pvectorTmp_trg.Eta();
                double phi_trg = pvectorTmp_trg.Phi();
                double pt_trg = pvectorTmp_trg.Pt();
                
                double effweight_trg = effhisto->GetBinContent(effhisto->FindBin(eta_trg,pt_trg));
                
                for(int nass=0;nass<nMult_ass;nass++)
                {
                    TVector3 pvectorTmp_ass = pVectTmp_ass[nass];
                    double eta_ass = pvectorTmp_ass.Eta();
                    double phi_ass = pvectorTmp_ass.Phi();
                    double pt_ass = pvectorTmp_ass.Pt();
                    
                    double effweight_ass = effhisto->GetBinContent(effhisto->FindBin(eta_ass,pt_ass));
                    
                    double deltaEta=eta_ass-eta_trg;
                    double deltaPhi=phi_ass-phi_trg;
                    if(deltaPhi>PI)
                        deltaPhi=deltaPhi-2*PI;
                    if(deltaPhi<-PI)
                        deltaPhi=deltaPhi+2*PI;
                    if(deltaPhi>-PI && deltaPhi<-PI/2.)
                        deltaPhi=deltaPhi+2*PI;
                    
                    //if(deltaEta==0 && deltaPhi==0) continue;
                    //if(fabs(deltaEta)<0.09 && fabs(deltaPhi)<0.09) continue;
                    
                    hBackground->Fill(deltaEta,deltaPhi,1.0/nMult_trg_eff/effweight_trg/effweight_ass);
                }
            }
        }
    }
    
    /*int nEvent = hMult->Integral(3,10000);
     double Bz = hBackground->GetBinContent(275);
     hCorrelation->Add(hSignal);
     hCorrelation->Scale(Bz);
     hCorrelation->Divide(hBackground);
     hCorrelation->Scale(1.0/nEvent);*/
}

//define this as a plug-in
DEFINE_FWK_MODULE(CorrAnaExampleMB);
