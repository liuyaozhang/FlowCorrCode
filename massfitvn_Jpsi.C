#include <iostream>

#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TString.h"
#include "TMath.h"

#include <vector>
#include <cmath>

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"

int iparmassfit_poly3bkg_floatwidth[11] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
int iparvnfit_poly3bkg_floatwidth[15] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};

struct GlobalChi2_poly3bkg_floatwidth {
    GlobalChi2_poly3bkg_floatwidth(  ROOT::Math::IMultiGenFunction & f1,
                                   ROOT::Math::IMultiGenFunction & f2) :
    fChi2_1(&f1), fChi2_2(&f2) {}
    
    // parameter vector is first background (in common 1 and 2)
    // and then is signal (only in 2)
    double operator() (const double *par) const {
        double p1[11];
        for(int i = 0; i < 11; ++i) p1[i] = par[iparmassfit_poly3bkg_floatwidth[i]];
        
        double p2[15];
        for(int i = 0; i < 15; ++i) p2[i] = par[iparvnfit_poly3bkg_floatwidth[i]];
        
        return (*fChi2_1)(p1) + (*fChi2_2)(p2);
    }
    
    const  ROOT::Math::IMultiGenFunction * fChi2_1;
    const  ROOT::Math::IMultiGenFunction * fChi2_2;
};

double crystalball_function(double x, double alpha, double n, double sigma, double mean) {
    // evaluate the crystal ball function
    if (sigma < 0.)     return 0.;
    double z = (x - mean)/sigma;
    if (alpha < 0) z = -z;
    double abs_alpha = std::abs(alpha);
    // double C = n/abs_alpha * 1./(n-1.) * std::exp(-alpha*alpha/2.);
    // double D = std::sqrt(M_PI/2.)*(1.+ROOT::Math::erf(abs_alpha/std::sqrt(2.)));
    // double N = 1./(sigma*(C+D));
    if (z  > - abs_alpha)
        return std::exp(- 0.5 * z * z);
    else {
        //double A = std::pow(n/abs_alpha,n) * std::exp(-0.5*abs_alpha*abs_alpha);
        double nDivAlpha = n/abs_alpha;
        double AA =  std::exp(-0.5*abs_alpha*abs_alpha);
        double B = nDivAlpha -abs_alpha;
        double arg = nDivAlpha/(B-z);
        return AA * std::pow(arg,n);
    }
}

/*p definitions
 [0] CB1 yield;
 [1] Common mean of CB and Gaus;
 [2] CB1 sigma;
 [3] CB n;
 [4] CB alpha;
 [5] CB2 yield;
 [6] CB2 sigma;
 [7-10] poly 3;
 [11] v2 signal;
 [12-13] v2 bkg;
 */
double crystalball_function_total(const double *x, const double *p) {
    // if ((!x) || (!p)) return 0.; // just a precaution
    // [Constant] * ROOT::Math::crystalball_function(x, [Alpha], [N], [Sigma], [Mean])
    return (p[0] * crystalball_function(x[0], p[3], p[4], p[2], p[1]) + p[5]*crystalball_function(x[0], p[3], p[4], p[6], p[1]) + p[7] + p[8]*x[0] + p[9]*x[0]*x[0] + p[10]*x[0]*x[0]*x[0]);
}

double crystalball_function_signal(const double *x, const double *p) {
    // if ((!x) || (!p)) return 0.; // just a precaution
    // [Constant] * ROOT::Math::crystalball_function(x, [Alpha], [N], [Sigma], [Mean])
    return (p[0] * crystalball_function(x[0], p[3], p[4], p[2], p[1]) + p[5]*crystalball_function(x[0], p[3], p[4], p[6], p[1]));
}

double crystalball_function_v2(const double *x, const double *p) {
    // if ((!x) || (!p)) return 0.; // just a precaution
    // [Constant] * ROOT::Math::crystalball_function(x, [Alpha], [N], [Sigma], [Mean])
    return ( ( p[0] * crystalball_function(x[0], p[3], p[4], p[2], p[1]) + p[5]*crystalball_function(x[0], p[3], p[4], p[6], p[1]) )/(p[0] * crystalball_function(x[0], p[3], p[4], p[2], p[1]) + p[5]*crystalball_function(x[0], p[3], p[4], p[6], p[1]) + p[7] + p[8]*x[0] + p[9]*x[0]*x[0] + p[10]*x[0]*x[0]*x[0]) * p[11] + (1 - ( p[0] * crystalball_function(x[0], p[3], p[4], p[2], p[1]) + p[5]*crystalball_function(x[0], p[3], p[4], p[6], p[1]) )/(p[0] * crystalball_function(x[0], p[3], p[4], p[2], p[1]) + p[5]*crystalball_function(x[0], p[3], p[4], p[6], p[1]) + p[7] + p[8]*x[0] + p[9]*x[0]*x[0] + p[10]*x[0]*x[0]*x[0])) * (p[12] + (exp(p[13] + p[14]*x[0]))) );
}

double function_v2_sig(const double *x, const double *p){
    return ( p[0] * crystalball_function(x[0], p[3], p[4], p[2], p[1]) + p[5]*crystalball_function(x[0], p[3], p[4], p[6], p[1]) )/(p[0] * crystalball_function(x[0], p[3], p[4], p[2], p[1]) + p[5]*crystalball_function(x[0], p[3], p[4], p[6], p[1]) + p[7] + p[8]*x[0] + p[9]*x[0]*x[0] + p[10]*x[0]*x[0]*x[0]) * p[11];
}

void massfitvn_Jpsi()
{
    double fit_range_low = 2.6;
    double fit_range_high = 3.5;
    double JPsi_mass = 3.097;
    int npt = 7;
    TFile* file1 = TFile::Open("HM185_JpsivnHist_etagap1p5_v30_eff_extdeta.root");
    
    TFile ofile("v2vspt_fromfit_jpsi_HM185_250_deta1p5_doubleCB_v30_eff_exp_extdeta.root","RECREATE");
    
    //v12
    double alpha_fit[14] = {4.30986,3.50841,3.03436,2.73741,2.37934,2.10685,2.03615};
    double n_fit[14] = {1.88853,1.9839,2.03198,2.07295,2.11001,2.15234,2.10154};
    
    TF1* fmasssig[9];
    TF1* fmassbkg[9];
    TF1* fmasstotal[9];
    TF1* fvn[9];
    
    double pt[13];
    double KET_ncq[13];
    double v2[13];
    double v2e[13];
    double v2_bkg[13];
    double v2_ncq[13];
    double v2e_ncq[13];
    double ptbin[14] = {0.2, 1.8, 3.0, 4.5, 6.0, 8.0, 10, 20};
    double a[13];
    double b[13];
    double sigfrac[13];
    
    TCanvas* c[10];
    for(int i=0;i<npt;i++)
    {
        c[i] = new TCanvas(Form("c_%d",i),Form("c_%d",i),800,400);
        c[i]->Divide(2,1);
    }
    
    for(int i=0;i<npt;i++)
    {
        c[i]->cd(1)->SetTopMargin(0.06);
        c[i]->cd(1)->SetLeftMargin(0.18);
        c[i]->cd(1)->SetRightMargin(0.043);
        c[i]->cd(1)->SetBottomMargin(0.145);
        c[i]->cd(2)->SetTopMargin(0.06);
        c[i]->cd(2)->SetLeftMargin(0.18);
        c[i]->cd(2)->SetRightMargin(0.043);
        c[i]->cd(2)->SetBottomMargin(0.145);

    }
    
    TCanvas* c2 = new TCanvas("c2","c2",100,100);
    
    TLatex* tex = new TLatex;
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.045);
    tex->SetLineWidth(2);
 
    TLatex* texCMS = new TLatex;
    texCMS->SetNDC();
    texCMS->SetTextFont(42);
    texCMS->SetTextSize(0.05);
    texCMS->SetTextAlign(12);
    
    TH1D* hist = new TH1D("hist","",10,2.6,3.5);
    hist->SetLineWidth(0);
    //hist->GetYaxis()->SetRangeUser(0,0.3);
    hist->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV)");
    hist->GetYaxis()->SetTitle("v_{2}^{S+B}");
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
    hist->GetXaxis()->SetTitleOffset(1.3);
    hist->GetYaxis()->SetTitleOffset(2);
    hist->GetXaxis()->SetLabelOffset(0.007);
    hist->GetYaxis()->SetLabelOffset(0.007);
    hist->GetXaxis()->SetTitleSize(0.045);
    hist->GetYaxis()->SetTitleSize(0.045);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetLabelSize(0.04);
    hist->SetMinimum(0.01);
    hist->SetMaximum(0.33);
    
    c2->cd();
    hist->Draw();
    
    for(int i=0;i<npt;i++)
    {
        TH1D* h_data = (TH1D*)file1->Get(Form("massjpsi_pt%d",i));
        h_data->SetMinimum(0);
        h_data->SetMarkerSize(0.8);
        h_data->SetMarkerStyle(20);
        h_data->SetLineWidth(1);
        h_data->SetOption("e");
        
        h_data->Rebin(2);

        h_data->GetXaxis()->SetRangeUser(2.6,3.5);
        h_data->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV)");
        h_data->GetYaxis()->SetTitle("Entries / 10 MeV");
        h_data->GetXaxis()->CenterTitle();
        h_data->GetYaxis()->CenterTitle();
        h_data->GetXaxis()->SetTitleOffset(1.3);
        h_data->GetYaxis()->SetTitleOffset(2);
        h_data->GetXaxis()->SetLabelOffset(0.007);
        h_data->GetYaxis()->SetLabelOffset(0.007);
        h_data->GetXaxis()->SetTitleSize(0.045);
        h_data->GetYaxis()->SetTitleSize(0.045);
        h_data->GetXaxis()->SetTitleFont(42);
        h_data->GetYaxis()->SetTitleFont(42);
        h_data->GetXaxis()->SetLabelFont(42);
        h_data->GetYaxis()->SetLabelFont(42);
        h_data->GetXaxis()->SetLabelSize(0.04);
        h_data->GetYaxis()->SetLabelSize(0.04);
        
        h_data->GetXaxis()->SetNoExponent(true);
        ((TGaxis*)h_data->GetXaxis())->SetMaxDigits(7);
        
        h_data->SetMaximum(h_data->GetMaximum()*1.5);
        
        TH1D* h_pt = (TH1D*)file1->Get(Form("Ptjpsi_eff_pt%d",i));
        TH1D* h_KET = (TH1D*)file1->Get(Form("KETjpsi_eff_pt%d",i));
        pt[i] = h_pt->GetMean();
        KET_ncq[i] = h_KET->GetMean()/2.0;

        c[i]->cd(1);
        
        /*p definitions
         [0] CB1 yield;
         [1] Common mean of CB and Gaus;
         [2] CB1 sigma;
         [3] CB n;
         [4] CB alpha;
         [5] CB2 yield;
         [6] CB2 sigma;
         [7-10] poly 3;
         [11] v2 signal;
         [12-13] v2 bkg;
         */
        TF1* f = new TF1(Form("f_%d",i), crystalball_function_total, fit_range_low, fit_range_high, 11);
        f->SetLineColor(2);
        f->SetLineWidth(1);
        f->SetParNames("CB1_Yield","common_mean","CB1_sigma","CB_N","CB_Alpha","CB2_Yield","CB2_Sigma","Pol0","Pol1","Pol2","Pol3");

        //first fit data mass signal + bkg
        
        f->SetParameter(0,10000.);
        f->SetParameter(1,JPsi_mass);
        f->SetParameter(2,0.03);
        f->SetParameter(3,1.0);
        f->SetParameter(4,1.0);
        f->SetParameter(5,10000);
        f->SetParameter(6,0.03);
        
        f->SetParLimits(2,0.01,0.1);
        f->SetParLimits(6,0.01,0.1);
        
        //fix alpha & n from MC
        f->FixParameter(4,alpha_fit[i]);
        f->FixParameter(3,n_fit[i]);
        
        f->FixParameter(1,JPsi_mass); //for first few attempt fix mean of gaussian to get reasonable estimation of other pars; later open it up
        h_data->Fit(Form("f_%d",i),"q","",fit_range_low,fit_range_high);
        h_data->Fit(Form("f_%d",i),"q","",fit_range_low,fit_range_high);
        f->ReleaseParameter(1); //now let gaussian mean float
        h_data->Fit(Form("f_%d",i),"L q","",fit_range_low,fit_range_high);
        h_data->Fit(Form("f_%d",i),"L q","",fit_range_low,fit_range_high);
        h_data->Fit(Form("f_%d",i),"L m","",fit_range_low,fit_range_high);
        
        
        //draw D0 signal separately
        TF1* f1 = new TF1(Form("f_sig_%d",i), crystalball_function_signal, fit_range_low, fit_range_high, 7);
        f1->SetLineColor(kOrange-3);
        f1->SetLineWidth(1);
        f1->SetLineStyle(2);
        f1->SetFillColorAlpha(kOrange-3,0.3);
        f1->SetFillStyle(1001);
        f1->FixParameter(0,f->GetParameter(0));
        f1->FixParameter(1,f->GetParameter(1));
        f1->FixParameter(2,f->GetParameter(2));
        f1->FixParameter(3,f->GetParameter(3));
        f1->FixParameter(4,f->GetParameter(4));
        f1->FixParameter(5,f->GetParameter(5));
        f1->FixParameter(6,f->GetParameter(6));
        
        fmasssig[i] = (TF1*)f1->Clone();
        fmasssig[i]->SetName(Form("masssigfcn_pt%d",i));
        fmasssig[i]->Write();
        
        f1->Draw("LSAME");
        
        //draw poly bkg separately
        TF1* f3 = new TF1(Form("f_bkg_%d",i),"[7] + [8]*x + [9]*x*x + [10]*x*x*x", fit_range_low, fit_range_high);
        f3->SetLineColor(4);
        f3->SetLineWidth(1);
        f3->SetLineStyle(2);
        f3->FixParameter(7,f->GetParameter(7));
        f3->FixParameter(8,f->GetParameter(8));
        f3->FixParameter(9,f->GetParameter(9));
        f3->FixParameter(10,f->GetParameter(10));
        
        fmassbkg[i] = (TF1*)f3->Clone();
        fmassbkg[i]->SetName(Form("massbkgfcn_pt%d",i));
        fmassbkg[i]->Write();
        
        f3->Draw("LSAME");
        
        tex->DrawLatex(0.22,0.86,"185 #leq N_{trk}^{offline} < 250");
        tex->DrawLatex(0.22,0.80,Form("%.1f < p_{T} < %.1f GeV",ptbin[i],ptbin[i+1]));
        tex->DrawLatex(0.22,0.74,"-2.86 < y_{cm} < -1.86 or 0.94 < y_{cm} < 1.94");
        
        texCMS->DrawLatex(.18,.97,"#font[61]{CMS} #it{Preliminary}");
        //texCMS->DrawLatex(.18,.97,"#font[61]{CMS}");
        texCMS->DrawLatex(0.73,0.97, "#scale[0.8]{pPb 8.16 TeV}");
        
        TLegend* leg = new TLegend(0.21,0.4,0.5,0.65,NULL,"brNDC");
        leg->SetBorderSize(0);
        leg->SetTextSize(0.045);
        leg->SetTextFont(42);
        leg->SetFillStyle(0);
        leg->AddEntry(h_data,"data","p");
        leg->AddEntry(f,"Fit","L");
        leg->AddEntry(f1,"J/#psi Signal","f");
        leg->AddEntry(f3,"Combinatorial","l");
        leg->Draw("SAME");
        
        sigfrac[i] = f1->Integral(2.94,3.24)/f->Integral(2.94,3.24);
        
        //c->Print(Form("plots/massfit_pt%d.pdf",i));
        
        //fit vn
        //[9] is vn_sig
        //[10-11] is vn bkg, const + linear vn(pT)
        TGraphErrors* vn_data = (TGraphErrors*)file1->Get(Form("v2_mass_pt%d",i));
        
        c[i]->cd(2);
        
        hist->Draw();
        
        TF1* fmass_combinemassvnfit = new TF1(Form("fmass_combinemassvnfit_%d",i),crystalball_function_total, fit_range_low, fit_range_high, 11);
        
        TF1* fvn_combinemassvnfit = new TF1(Form("fvn_combinemassvnfit_%d",i), crystalball_function_v2, fit_range_low, fit_range_high, 15);
        
        fmass_combinemassvnfit->SetLineColor(2);
        fmass_combinemassvnfit->SetLineWidth(1);
        
        fvn_combinemassvnfit->SetLineColor(2);
        fvn_combinemassvnfit->SetLineWidth(1);

        ROOT::Math::WrappedMultiTF1 wfmass_combinemassvnfit(*fmass_combinemassvnfit,1);
        ROOT::Math::WrappedMultiTF1 wfvn_combinemassvnfit(*fvn_combinemassvnfit,1);
        
        ROOT::Fit::DataOptions opt;
        ROOT::Fit::DataRange range_massfit;

        range_massfit.SetRange(fit_range_low,fit_range_high);
        ROOT::Fit::BinData datamass(opt,range_massfit);
        ROOT::Fit::FillData(datamass, h_data);
        
        ROOT::Fit::DataRange range_vnfit;
        range_vnfit.SetRange(fit_range_low,fit_range_high);
        ROOT::Fit::BinData datavn(opt,range_vnfit);
        ROOT::Fit::FillData(datavn, vn_data);
        
        ROOT::Fit::Chi2Function chi2_B(datamass, wfmass_combinemassvnfit);
        ROOT::Fit::Chi2Function chi2_SB(datavn, wfvn_combinemassvnfit);
        
        GlobalChi2_poly3bkg_floatwidth globalChi2(chi2_B, chi2_SB);

        ROOT::Fit::Fitter fitter;
        
        const int Npar = 15;
        double par0[Npar];
        for( int ipar = 0; ipar < f->GetNpar(); ipar++ ) par0[ipar] = f->GetParameter(ipar);
        par0[11] = 0.01;
        par0[12] = 0.10;
        par0[13] = 0.05;
        par0[14] = 0.01;

        
        fitter.Config().SetParamsSettings(Npar,par0);
        // fix parameter
        fitter.Config().ParSettings(0).Fix();
        fitter.Config().ParSettings(1).Fix();
        fitter.Config().ParSettings(2).Fix();
        fitter.Config().ParSettings(3).Fix();
        fitter.Config().ParSettings(4).Fix();
        fitter.Config().ParSettings(5).Fix();
        fitter.Config().ParSettings(6).Fix();
        fitter.Config().ParSettings(7).Fix();
        fitter.Config().ParSettings(8).Fix();
        fitter.Config().ParSettings(9).Fix();
        fitter.Config().ParSettings(10).Fix();
        
        fitter.Config().MinimizerOptions().SetPrintLevel(0);
        fitter.Config().SetMinimizer("Minuit2","Migrad");

        fitter.FitFCN(Npar,globalChi2,0,datamass.Size()+datavn.Size(),true);
        ROOT::Fit::FitResult result = fitter.Result();
        result.Print(std::cout);
        
        fmass_combinemassvnfit->SetFitResult( result, iparmassfit_poly3bkg_floatwidth);
        fmass_combinemassvnfit->SetRange(range_massfit().first, range_massfit().second);
        fmass_combinemassvnfit->SetLineColor(kRed);
        h_data->GetListOfFunctions()->Add(fmass_combinemassvnfit);
        //c->cd();
        //h_data->Draw();
        
        fvn_combinemassvnfit->SetFitResult( result, iparvnfit_poly3bkg_floatwidth);
        fvn_combinemassvnfit->SetRange(range_vnfit().first, range_vnfit().second);
        fvn_combinemassvnfit->SetLineColor(2);
        //fvn_combinemassvnfit->SetLineStyle(2);
        vn_data->GetListOfFunctions()->Add(fvn_combinemassvnfit);
        vn_data->SetTitle("");
        vn_data->SetMarkerSize(0.8);
        vn_data->SetLineWidth(1);
        //c1->cd();
        vn_data->Draw("PESAME");
        
        fvn[i] = (TF1*)fvn_combinemassvnfit->Clone();
        fvn[i]->SetName(Form("vnfit_pt%d",i));
        fvn[i]->Write();
        
        fmasstotal[i] = (TF1*)fmass_combinemassvnfit->Clone();
        fmasstotal[i]->SetName(Form("masstotalfcn_pt%d",i));
        fmasstotal[i]->Write();
        
        tex->DrawLatex(0.22,0.86,"185 #leq N_{trk}^{offline} < 250");
        tex->DrawLatex(0.22,0.80,Form("%.1f < p_{T} < %.1f GeV",ptbin[i],ptbin[i+1]));
        //tex->DrawLatex(0.22,0.74,"1.4 < |y_{cm}+0.46| < 2.4");
        tex->DrawLatex(0.22,0.74,"-2.86 < y_{cm} < -1.86 or 0.94 < y_{cm} < 1.94");
        //tex->DrawLatex(0.22,0.68,"|#Delta#eta| > 2");

        
        //texCMS->DrawLatex(.18,.97,"#font[61]{CMS}");
        texCMS->DrawLatex(.18,.97,"#font[61]{CMS} #it{Preliminary}");
        texCMS->DrawLatex(0.73,0.97, "#scale[0.8]{pPb 8.16 TeV}");
        
        v2[i] = fvn_combinemassvnfit->GetParameter(11);
        v2e[i] = fvn_combinemassvnfit->GetParError(11);
        v2_bkg[i] = fvn_combinemassvnfit->GetParameter(12) + fvn_combinemassvnfit->GetParameter(13) * JPsi_mass;
        v2_ncq[i] = v2[i]/2.0;
        v2e_ncq[i] = v2e[i]/2.0;
        a[i] = fvn_combinemassvnfit->GetParameter(12);
        b[i] = fvn_combinemassvnfit->GetParameter(13);
        
        TF1* fvnbkg = new TF1(Form("fvnbkg_%d",1),"( [0] + [1] * x)", fit_range_low, fit_range_high);
        fvnbkg->FixParameter(0,fvn_combinemassvnfit->GetParameter(12));
        fvnbkg->FixParameter(1,fvn_combinemassvnfit->GetParameter(13));
        
        fvnbkg->SetName(Form("fvnbkg_fcn_pt%d",i));
        fvnbkg->Write();
        
        fvnbkg->SetLineStyle(7);
        //fvnbkg->Draw("LSAME");
        
        TF1* fvnsig = new TF1(Form("fvnsig_%d",i),function_v2_sig,fit_range_low,fit_range_high,12);
        for(int k=0;k<12;k++)
        {
            fvnsig->FixParameter(k,fvn_combinemassvnfit->GetParameter(k));
            
        }
        
        fvnsig->SetLineColor(kOrange-3);
        fvnsig->SetLineWidth(1);
        fvnsig->SetLineStyle(2);
        fvnsig->SetFillColorAlpha(kOrange-3,0.3);
        fvnsig->SetFillStyle(1001);
        
        //fvnsig->Draw("LSAME");
        
        TLegend* leg1 = new TLegend(0.72,0.525,0.91,0.65,NULL,"brNDC");
        leg1->SetBorderSize(0);
        leg1->SetTextSize(0.045);
        leg1->SetTextFont(42);
        leg1->SetFillStyle(0);
        leg1->AddEntry(h_data,"data","p");
        leg1->AddEntry(fvn_combinemassvnfit,"Fit","l");
        //leg1->AddEntry(fvnsig,"#alpha(#it{m}_{#mu#mu})v_{2}^{S}","f");
        leg1->Draw("SAME");

        double xmass[200];
        double pullmass[200];
        
        float Chi2=0;
        int ndf = (fit_range_high - fit_range_low)/0.01 - 8;
        
        for(int k=0;k<h_data->GetNbinsX();k++)
        {
            xmass[k] = h_data->GetBinCenter(k);
            pullmass[k] = (h_data->GetBinContent(k) - fmass_combinemassvnfit->Eval(xmass[k]))/h_data->GetBinError(k);
            if(fabs(pullmass[k])<5)
            {
                //cout<<pullmass[k]<<endl;
                Chi2 += pullmass[k]*pullmass[k];
            }
        }

        c[i]->cd(1);
        tex->DrawLatex(0.22,0.67,Form("#chi^{2}/ndf = %.0f/%d",Chi2,ndf));
        
        double xv2[200];
        double pullv2[200];
        double v2y[200];
        
        float Chi2v2=0;
        int ndfv2 = 8 - 4; //Nbin - Npar
        
        for(int k=0;k<vn_data->GetN()-1;k++)
        {
            vn_data->GetPoint(k,xv2[k],v2y[k]);
            //xv2[k] = vn_dara->GetBinCenter(k);
            pullv2[k] = (v2y[k] - fvn_combinemassvnfit->Eval(xv2[k]))/vn_data->GetErrorY(k);
            cout<<k<<": "<<pullv2[k]<<endl;
            if(fabs(pullv2[k])<1000)
            {
                //cout<<pullmass[k]<<endl;
                Chi2v2 += pullv2[k]*pullv2[k];
            }
            cout<<"fcn: "<<fvn_combinemassvnfit->Eval(xv2[k])<<endl;
            cout<<"data: "<<v2y[k]<<endl;
        }

        c[i]->cd(2);
        tex->DrawLatex(0.22,0.67,Form("#chi^{2}/ndf = %.1f/%d",Chi2v2,ndfv2));
        
    }
    
    for(int i=0;i<npt;i++)
    {
        c[i]->Print(Form("plots/v30/eff/exp/JPsi_mass_vnfit_combine_pt%d.pdf",i));
        c[i]->Print(Form("plots/v30/eff/exp/JPsi_mass_vnfit_combine_pt%d.gif",i));
    }
    
    TGraphErrors* v2plot = new TGraphErrors(npt,pt,v2,0,v2e);
    TGraphErrors* v2ncqplot = new TGraphErrors(npt,KET_ncq,v2_ncq,0,v2e_ncq);
    TGraphErrors* v2bkgplot = new TGraphErrors(npt,pt,v2_bkg,0,0);
    
    v2plot->SetName("v2vspt");
    v2ncqplot->SetName("v2vsKET_ncq");
    v2bkgplot->SetName("v2bkgvspt");
    
    v2plot->Write();
    v2ncqplot->Write();
    v2bkgplot->Write();
}
