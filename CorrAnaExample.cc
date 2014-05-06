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
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>
//
// class declaration
//

#define PI 3.1416

#define effkss(x) ((1.09925e+11*TMath::Power(x,1.50511e+01)*TMath::Exp(6.58074*x*x-2.94487e+01*x+9.72194))/(6.71504e+11*TMath::Power(x,3.19081)*TMath::Exp(1.97717*x*x-9.90447*x-4.65781)))

#define effksb(x) ((6.47559e+02*TMath::Power(x,2.95369e-01)*TMath::Exp(9.18237e-02*x*x-1.89678*x+7.63891))/(1.20601e+01*TMath::Power(x,-1.86165)*TMath::Exp(4.17033e-02*x*x-9.39659e-01*x+1.30386e+01)))

#define efflas(x) ((4.76443e+06*TMath::Power(x,1.62889e+01)*TMath::Exp(2.74004*x*x-1.97581e+01*x+1.16309e+01))/(5.30771e+08*TMath::Power(x,3.59273)*TMath::Exp(1.50703*x*x-8.45701*x+9.43797e-01)))

#define efflab(x) ((3.86297e-01*TMath::Power(x,1.91207)*TMath::Exp(8.37588e-02*x*x-2.15583*x+1.33689e+01))/(7.01220*TMath::Power(x,-4.80662e-01)*TMath::Exp(7.33837e-02*x*x-1.53854*x+1.35978e+01)))


using namespace std;

class CorrAnaExample : public edm::EDAnalyzer {
public:
    explicit CorrAnaExample(const edm::ParameterSet&);
    ~CorrAnaExample();
    
    
private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;
    
    // ----------member data ---------------------------
    TH1D* hMult_selected;
    TH1D* hMult;
    TH2D* hSignal;
    TH2D* hBackground;
    //TH2D* hCorrelation;
    
    TH1D* hMass_ks[13];
    TH1D* hMass_la[13];
    TH1D* hMass_Xi[13];
    TH1D* hMass_Omg[13];
    
    TH2F* effhisto;
    
    double ptcut[14];
    
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
    edm::Service<TFileService> fs;
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
CorrAnaExample::CorrAnaExample(const edm::ParameterSet& iConfig)

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
}


CorrAnaExample::~CorrAnaExample()
{
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
CorrAnaExample::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    
    double ptcut[14] = {0.2,0.4,0.6,0.8,1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,9.0,12.0};
    
    // select on requirement of valid vertex
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByLabel("offlinePrimaryVertices",vertices);
    double bestvz=-999.9, bestvx=-999.9, bestvy=-999.9;
    double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
    const reco::Vertex & vtx = (*vertices)[0];
    bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
    bestvzError = vtx.zError(); bestvxError = vtx.xError(); bestvyError = vtx.yError();
    
    if(bestvz < -15.0 || bestvz>15.0) return;
    
    edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates_ks;
    iEvent.getByLabel("generalV0CandidatesNew","Kshort",v0candidates_ks);
    if(!v0candidates_ks.isValid()) return;
    
    edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates_la;
    iEvent.getByLabel("generalV0CandidatesNew","Lambda",v0candidates_la);
    if(!v0candidates_la.isValid()) return;
    
    /*edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates_xi;
    iEvent.getByLabel("generalV0CandidatesNew","Xi",v0candidates_xi);
    if(!v0candidates_xi.isValid()) return;
    
    edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates_omg;
    iEvent.getByLabel("generalV0CandidatesNew","Omega",v0candidates_omg);
    if(!v0candidates_omg.isValid()) return;*/
    
    pVect_trg = new vector<TVector3>;
    pVect_ass = new vector<TVector3>;
    //----- loop over tracks -----
    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByLabel("generalTracks",tracks);
    
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
    hMult_selected->Fill(nMult_ass_good);
   
    if(nMult_ass_good<multMax_ && nMult_ass_good>=multMin_){
        for(unsigned it=0; it<v0candidates_ks->size(); ++it){
            
            const reco::VertexCompositeCandidate & trk = (*v0candidates_ks)[it];
            
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
            
            TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
            TVector3 secvec(px,py,pz);
            
            double agl = cos(secvec.Angle(ptosvec));
            
            if(agl<=0.999) continue;
            
            typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
            typedef ROOT::Math::SVector<double, 3> SVector3;
            
            SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
            SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);
            
            double dl = ROOT::Math::Mag(distanceVector);
            double dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl;
            
            double dlos = dl/dlerror;
            
            if(dlos<=5) continue;
            
            const reco::Candidate * dau1 = trk.daughter(0);
            const reco::Candidate * dau2 = trk.daughter(1);
            
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
            
            TVector3 dauvecsum(dauvec1+dauvec2);
            double v0masspiproton1 = sqrt((sqrt(0.93827*0.93827+pd1*pd1)+sqrt(0.13957*0.13957+pd2*pd2))*(sqrt(0.93827*0.93827+pd1*pd1)+sqrt(0.13957*0.13957+pd2*pd2))-dauvecsum.Mag2());
            
            double v0masspiproton2 = sqrt((sqrt(0.13957*0.13957+pd1*pd1)+sqrt(0.93827*0.93827+pd2*pd2))*(sqrt(0.13957*0.13957+pd1*pd1)+sqrt(0.93827*0.93827+pd2*pd2))-dauvecsum.Mag2());
            
            if((v0masspiproton1>=(1.115683-0.010) && v0masspiproton1<=(1.115683+0.010)) || (v0masspiproton2>=(1.115683-0.010) && v0masspiproton2<=(1.115683+0.010)) ) continue;
            
            //efficiency
            double effks = 1.0;
            
            for(int i=0;i<13;i++)
            {
                if(eta<=etaMax_trg_ && eta>=etaMin_trg_ && pt<=ptcut[i+1] && pt>=ptcut[i]){
                    hMass_ks[i]->Fill(mass,1.0/effks);
                }
            }
        }
        
        for(unsigned it=0; it<v0candidates_la->size(); ++it){
            
            const reco::VertexCompositeCandidate & trk = (*v0candidates_la)[it];
            
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
            
            TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
            TVector3 secvec(px,py,pz);
            
            double agl = cos(secvec.Angle(ptosvec));
            
            if(agl<=0.999) continue;
            
            typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
            typedef ROOT::Math::SVector<double, 3> SVector3;
            
            SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
            SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);
            
            double dl = ROOT::Math::Mag(distanceVector);
            double dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl;
            
            double dlos = dl/dlerror;
            
            if(dlos<=5) continue;
            
            const reco::Candidate * dau1 = trk.daughter(0);
            const reco::Candidate * dau2 = trk.daughter(1);
            
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
            
            TVector3 dauvecsum(dauvec1+dauvec2);
            double v0masspipi = sqrt((sqrt(0.13957*0.13957+pd1*pd1)+sqrt(0.13957*0.13957+pd2*pd2))*(sqrt(0.13957*0.13957+pd1*pd1)+sqrt(0.13957*0.13957+pd2*pd2))-dauvecsum.Mag2());
            double v0massee = sqrt((sqrt(0.000511*0.000511+pd1*pd1)+sqrt(0.000511*0.000511+pd2*pd2))*(sqrt(0.000511*0.000511+pd1*pd1)+sqrt(0.000511*0.000511+pd2*pd2))-dauvecsum.Mag2());
            
            if( (v0masspipi>=(0.497614-0.020) && v0masspipi<=(0.497614+0.020)) || v0massee <= 0.015 ) continue;
            
            //efficiency
            double effla = 1.0;
            
            for(int i=0;i<13;i++)
            {
                if(eta<=etaMax_trg_ && eta>=etaMin_trg_ && pt<=ptcut[i+1] && pt>=ptcut[i]){
                    hMass_la[i]->Fill(mass,1.0/effla);
                }
            }
        }
        
        /*for(unsigned it=0; it<v0candidates_xi->size(); ++it){
            
            const reco::VertexCompositeCandidate & trk = (*v0candidates_xi)[it];
            
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
            
            TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
            TVector3 secvec(px,py,pz);
            
            double agl = cos(secvec.Angle(ptosvec));
            
            if(agl<=0.9999) continue;
            
            typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
            typedef ROOT::Math::SVector<double, 3> SVector3;
            
            SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
            SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);
            
            double dl = ROOT::Math::Mag(distanceVector);
            double dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl;
            
            double dlos = dl/dlerror;
            
            if(dlos<=0) continue;
            
            //efficiency
            double effxi = 1.0;
            
            for(int i=0;i<13;i++)
            {
                if(trk.eta()<=etaMax_trg_ && trk.eta()>=etaMin_trg_ && trk.pt()<=ptcut[i+1] && trk.pt()>=ptcut[i]){
                    hMass_Xi[i]->Fill(mass,1.0/effxi);
                }
            }
        }

        for(unsigned it=0; it<v0candidates_omg->size(); ++it){
            
            const reco::VertexCompositeCandidate & trk = (*v0candidates_omg)[it];
            
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
            
            TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
            TVector3 secvec(px,py,pz);
            
            double agl = cos(secvec.Angle(ptosvec));
            
            if(agl<=0.9999) continue;
            
            typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
            typedef ROOT::Math::SVector<double, 3> SVector3;
            
            SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
            SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);
            
            double dl = ROOT::Math::Mag(distanceVector);
            double dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl;
            
            double dlos = dl/dlerror;
            
            if(dlos<=0) continue;
            
            //efficiency
            double effomg = 1.0;
            
            for(int i=0;i<13;i++)
            {
                if(trk.eta()<=etaMax_trg_ && trk.eta()>=etaMin_trg_ && trk.pt()<=ptcut[i+1] && trk.pt()>=ptcut[i]){
                    hMass_Omg[i]->Fill(mass,1.0/effomg);
                }
            }
        }*/

        
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
        if(eta<=etaMax_trg_ && eta>=etaMin_trg_ && pt<=ptMax_trg_ && pt>=ptMin_trg_) pVect_trg->push_back(pvector);
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
            if(fabs(deltaEta)<0.028 && fabs(deltaPhi)<0.02) continue;
            hSignal->Fill(deltaEta,deltaPhi,1.0/nMult_trg_eff/effweight_trg/effweight_ass);
        }
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
CorrAnaExample::beginJob(){
    
    TH1D::SetDefaultSumw2();
    
    edm::FileInPath fip("Demo/DemoAnalyzer/data/TrackCorrections_HIJING_538_OFFICIAL_Mar24.root");
    TFile f(fip.fullPath().c_str(),"READ");
    effhisto = (TH2F*)f.Get("rTotalEff3D");

    hMult_selected = fs->make<TH1D>("mult_selected",";N",300,0,300);
    hMult = fs->make<TH1D>("mult",";N",1500,0,1500);
    hSignal = fs->make<TH2D>("signal",";#Delta#eta;#Delta#phi",
                             33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
    hBackground = fs->make<TH2D>("background",";#Delta#eta;#Delta#phi",
                                 33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
    //hCorrelation = fs->make<TH2D>("correlation",";#Delta#eta;#Delta#phi",
    //                              27,-4.05,4.05,33,-(0.5+1/32)*PI,(1.5+1/32)*PI);
    for(int i=0; i<13; i++)
    {
        hMass_ks[i] = fs->make<TH1D>(Form("masskshort_pt%d",i),";GeV",2000,0,1.0);
        hMass_la[i] = fs->make<TH1D>(Form("masslambda_pt%d",i),";GeV",2000,0.5,1.5);
        //hMass_Xi[i] = fs->make<TH1D>(Form("massXi_pt%d",i),";GeV",2000,0.7,1.7);
        //hMass_Omg[i] = fs->make<TH1D>(Form("massOmega_pt%d",i),";GeV",2000,1.0,2.0);
    }
    pVectVect_trg = new vector< vector<TVector3> >;
    pVectVect_ass = new vector< vector<TVector3> >;
    zvtxVect = new vector<double>;
}

// ------------ method called once each job just after ending the event loop  ------------
void
CorrAnaExample::endJob() {
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
                    if(fabs(deltaEta)<0.028 && fabs(deltaPhi)<0.02) continue;
                    
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
DEFINE_FWK_MODULE(CorrAnaExample);
