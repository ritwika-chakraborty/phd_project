

#include <TCanvas.h>

#include <TH1D.h>

#include <TF1.h>

#include <TMath.h>

#include <TRandom3.h>

#include <TGraphErrors.h>

#include <TFitResultPtr.h>

#include <TFitResult.h>

​

class BinIntegralFunc {

​

public:

  BinIntegralFunc(TF1* inputFunc, double inputBinWidth)

    : f(inputFunc)

    , binWidth(inputBinWidth)

  {}

​

  double Evaluate(double *x, double *p) {

    f->SetParameters(p);

    return f->Integral( x[0] - 0.5*binWidth, x[0] + 0.5*binWidth) / binWidth;

  }

​

private:

  TF1* f;

  double binWidth;

};

​

void tratiofitsim(){

​

  // Params

  int nEntries = 1000000000;

  int nTrials = 1;

  double maxTime = 450000;

  double binWidth = 150;

  int nBins = maxTime/binWidth;

​

  // Truth parameters

  double tau = 64000;

  double A = 0.42;

  double w_a = 0.00144;

  double phi = TMath::Pi();

​

  // Get N0 from integral of function (need to do numeric as we truncate exponential slightly and wiggle introduces error too)

  auto normalisedFunc = new TF1("truthFunc","[0]*TMath::Exp(-x/[1])*(1+[2]*cos([3]*x+[4]))",0,maxTime);

  normalisedFunc->SetParameters(1, tau, A, w_a, phi);

  double N0 = (binWidth*nEntries)/normalisedFunc->Integral(0,maxTime);

  

  // Declare truth function

  auto truthFunc = new TF1("truthFunc","[0]*TMath::Exp(-x/[1])*(1+[2]*cos([3]*x+[4]))",0,maxTime);

  truthFunc->SetParameters(N0, tau, A, w_a, phi);

  truthFunc->SetLineColor(2);

  truthFunc->SetLineWidth(1);

  truthFunc->SetNpx(10000);

​

  // Histograms to store pulls

  auto hA   = new TH1D("hA",";A Pull;Entries",500,-5,5);

  auto hw_a = new TH1D("hw_a",";#omega_{a} Pull;Entries",500,-5,5);

  auto hphi = new TH1D("hphi",";#phi Pull;Entries",500,-5,5);

  auto hConst = new TH1D("hConst",";Constant Term Pull;Entries",500,-5,5);

  auto hPVal = new TH1D("hPVal",";P-Value;Entries",100,0,1);

​

  // Do some pseudo data trials

  TRandom3* rng = new TRandom3(54321);

  double half_t_a = TMath::Pi()/w_a;

  double delta = 0; // Constant term

​

  // Normalisation for weighted parts

  double denom = 2+2*cosh(half_t_a/tau);

  double VThresh = 2./denom;

  double UThresh_p = VThresh + exp(+half_t_a/tau)/denom;

  double UThresh_m = UThresh_p + exp(-half_t_a/tau)/denom;

  std::cout << "Cut offs : " << VThresh << "\t" << UThresh_p << "\t" << UThresh_m << endl;

​

  for(int trial = 0; trial < nTrials; trial++){

    

    if(trial % 1 == 0) cout << "Trial " << trial << "/" << nTrials << endl;

​

    // Psuedo data histograms

    auto hData = new TH1F(Form("hData_%d",trial),";Time [ns]; Entries",nBins,0,maxTime);

    auto hData_U = new TH1F(Form("hData_U_%d",trial),";Time [ns]; Entries",nBins,0,maxTime);

    auto hData_V = new TH1F(Form("hData_V_%d",trial),";Time [ns]; Entries",nBins,0,maxTime);

    auto hData_num = new TH1F(Form("hData_num_%d",trial),";Time [ns]; Entries",nBins,0,maxTime);

    auto hData_den = new TH1F(Form("hData_den_%d",trial),";Time [ns]; Entries",nBins,0,maxTime);

​

    int nRandEntries = rng->Poisson(nEntries);

    for(int i = 0; i < nRandEntries; i++){

      double time = truthFunc->GetRandom();

      hData->Fill(time);

      double r = rng->Uniform();

      if(r < VThresh){

	hData_V->Fill(time);

      } else if (r < UThresh_p){

	hData_U->Fill(time+half_t_a);

      } else {

	hData_U->Fill(time-half_t_a);

      }

    }

​

    // Make ratio fit graph

    hData_num->Add(hData_U,hData_V,-1);

    hData_den->Add(hData_U,hData_V);

    auto gData_ratio = new TGraphErrors();

    for(int i = 1; i <= hData_num->GetNbinsX(); i++){

      int nPt = gData_ratio->GetN();

      if(hData_den->GetBinContent(i) > 0 && hData_num->GetBinCenter(i) > half_t_a){

	double ratio = hData_num->GetBinContent(i)/hData_den->GetBinContent(i);

	double ratioErr = sqrt((1. - ratio*ratio) / hData_den->GetBinContent(i));

	gData_ratio->SetPoint(nPt,hData_num->GetBinCenter(i),ratio);

	gData_ratio->SetPointError(nPt,0,ratioErr);

      }

    }

    

    // Fit function

    auto fitFunc = new TF1(Form("fitFunc_%d",trial),"[0]*cos([1]*x+[2])+[3]",half_t_a,100000);

    fitFunc->SetParameter(0, A);

    fitFunc->SetParameter(1, w_a);

    fitFunc->SetParameter(2, phi);

    fitFunc->SetParameter(3, delta);

    fitFunc->SetNpx(10000);

​

    // Fit function with integral over bin width

    auto binIntegralFunc = new BinIntegralFunc(fitFunc,binWidth);

    auto fitFuncInt = new TF1("f", binIntegralFunc, &BinIntegralFunc::Evaluate, fitFunc->GetXmin(), fitFunc->GetXmax(), fitFunc->GetNpar()); 

    fitFuncInt->SetParameters(fitFunc->GetParameters());

    fitFuncInt->SetParameter(3, fitFuncInt->GetParameter(3));

    fitFuncInt->SetNpx(10000);

    TFitResultPtr fitResult = gData_ratio->Fit(fitFuncInt,"RNSQ");

    fitResult->Print();

    

    // Get results from integral fit and set fitFunc parameter errors based on these

    fitFunc->SetParErrors(fitFuncInt->GetParErrors());

​

    // Fill histograms

    hA  ->Fill((fitFunc->GetParameter(0) - truthFunc->GetParameter(2))/fitFunc->GetParError(0));

    hw_a->Fill((fitFunc->GetParameter(1) - truthFunc->GetParameter(3))/fitFunc->GetParError(1));

    hphi->Fill((fitFunc->GetParameter(2) - truthFunc->GetParameter(4))/fitFunc->GetParError(2));

    hConst->Fill((fitFunc->GetParameter(3) - delta)/fitFunc->GetParError(3));

    hPVal->Fill(TMath::Prob(fitResult->Chi2(),fitResult->Ndf()));

    

    // Draw first one

    if(trial == 0){

      auto c1 = new TCanvas("c1","c1",200,10,800,600);

      hData->Draw();

      truthFunc->Draw("SAME");

      c1->SaveAs("Images/StartingData.png");

​

      auto c_U = new TCanvas("c_U","c_U",200,10,800,600);

      hData_U->Draw();

      c_U->SaveAs("Images/U.png");

​

      auto c_V = new TCanvas("c_V","c_V",200,10,800,600);

      hData_V->Draw();

      c_V->SaveAs("Images/V.png");

​

      auto c_num = new TCanvas("c_num","c_num",200,10,800,600);

      hData_num->Draw();

      c_num->SaveAs("Images/Numerator.png");

​

      auto c_den = new TCanvas("c_den","c_den",200,10,800,600);

      hData_den->Draw();

      c_den->SaveAs("Images/Denominator.png");

​

      auto c6 = new TCanvas("c6","c6",200,10,800,600);

      gData_ratio->Draw("AP");

      fitFuncInt->SetLineColor(4);

      fitFuncInt->SetLineWidth(1);

      fitFuncInt->Draw("SAME");

      fitFunc->SetLineColor(2);

      fitFunc->SetLineWidth(1);

      fitFunc->Draw("SAME");

      c6->SaveAs("Images/Ratio.png");

​

    } else {

      delete hData;

      delete hData_U;

      delete hData_V;

      delete hData_num;

      delete hData_den;

      delete gData_ratio;

      delete fitFunc;

      delete fitFuncInt;

      delete binIntegralFunc;

    }

​

  }

​

  auto cA = new TCanvas("cA","cA",200,10,800,600);

  hA->Draw();

  cA->SaveAs("Images/Pull_A.png");

​

  auto cw_a = new TCanvas("cw_a","cw_a",200,10,800,600);

  hw_a->Draw();

  cw_a->SaveAs("Images/Pull_wA.png");

​

  auto cphi = new TCanvas("cphi","cphi",200,10,800,600);

  hphi->Draw();

  cphi->SaveAs("Images/Pull_Phi.png");

​

  auto cConst = new TCanvas("cConst","cConst",200,10,800,600);

  hConst->Draw();

  cConst->SaveAs("Images/Pull_Const.png");

​

  auto cPVal = new TCanvas("cPVal","cPVal",200,10,800,600);

  hPVal->Draw();    

  cPVal->SaveAs("Images/PValue.png");

}

