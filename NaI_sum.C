#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TBranch.h>
#include <TH1I.h>
#include <TObjArray.h>
#include <TString.h>
#include <iostream>
#include <algorithm>

void ReadDataFile(TString filename = "./co_test_3_9_2024_0cm/UNFILTERED/SData_co_test_3_9_2024_0cm") {
  TFile *f = new TFile(filename + ".root", "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  } else {
    std::cout << "Reading file: " << filename << ".root" << std::endl;
  }

  TFile *fnew = new TFile(filename + "_processed.root", "RECREATE");

  Float_t hRange = 4096;
  const Int_t hBins = 1024;

  UShort_t Channel, Board, Energy;
  ULong64_t Timestamp;
  
  const Float_t gainUp = 6.05;
  const Float_t offsetUp = 12.51;
  const Float_t gainDown= 4.47;
  const Float_t offsetDown = 28.99;
  UShort_t CalibratedEnergy = 0.0;

  TNtuple *ntuple = (TNtuple *)f->Get("Data");
  if (!ntuple) {
    std::cerr << "TNtuple 'Data' not found in file." << std::endl;
    f->Close();
    return;
  }
  
  TNtuple *ntupleCal = new TNtuple("CalibratedData", "Calibrated Data", "Channel:Board:Timestamp:CalibratedEnergy");

  ntuple->SetBranchAddress("Channel", &Channel);
  ntuple->SetBranchAddress("Board", &Board);
  ntuple->SetBranchAddress("Timestamp", &Timestamp);
  ntuple->SetBranchAddress("Energy", &Energy);

  TH1I chUp("chUp", "Half Up spectrum", hBins, 0., hRange);
  chUp.GetXaxis()->SetTitle("Channels");

  TH1I chDown("chDown", "Half Down spectrum", hBins, 0., hRange);
  chDown.GetXaxis()->SetTitle("Channels");

  TH1I chVirtual("chVirtual", "Virtual Channel spectrum", hBins, 0., hRange);
  chVirtual.GetXaxis()->SetTitle("Channels");
  
  TH1F chUpCal("chUpCal", "Half Up spectrum (Calibrated)", hBins, 0., hRange);
  chUpCal.GetXaxis()->SetTitle("Energy");

  TH1F chDownCal("chDownCal", "Half Down spectrum (Calibrated)", hBins, 0., hRange);
  chDownCal.GetXaxis()->SetTitle("Energy");

  Long64_t nentries = ntuple->GetEntries();
  for (Long64_t i = 0; i < nentries; i++) {
    ntuple->GetEntry(i);

    if (Board <= 700) {
      if (Channel == 0) {
      	CalibratedEnergy = gainUp * Energy + offsetUp;
	chUpCal.Fill(CalibratedEnergy);
        chUp.Fill(Energy);
      } else if (Channel == 4) {
      	CalibratedEnergy = gainDown * Energy + offsetDown;
	chDownCal.Fill(CalibratedEnergy);
        chDown.Fill(Energy);
      }
    } else {
      chVirtual.Fill(Energy);
    }
    ntupleCal->Fill(Channel, Board, Timestamp, CalibratedEnergy);
  }

  std::cout << "Histograms created for Up Half, Down Half, and Virtual Channel." << std::endl;

  fnew->cd();
  chUp.Write();
  chDown.Write();
  chUpCal.Write();
  chDownCal.Write();
  chVirtual.Write();

  TNtuple *ntupleCopy = (TNtuple*)ntuple->CloneTree(-1); // Copy all entries
  ntupleCopy->SetName("Data");
  ntupleCopy->Write();
  std::cout << "TNtuple 'Data' copied to the new file." << std::endl;
  ntupleCal->Write();

  f->Close();
  fnew->Close();
}

//--------------------------------------------------------------------------------------------//

void CreateDTHistograms(TString filename = "./co_test_3_9_2024_0cm/UNFILTERED/SData_co_test_3_9_2024_0cm") {
    TFile *f = new TFile(filename + "_processed.root", "UPDATE");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    Float_t hRange = 20000;
    Double_t dTrange = 1e7;  
    Long64_t dTbins = dTrange / 1000;
    const Int_t hBins = 1024;

    TH1I channelsdt("channelsdt", "Segments dT spectrum", dTbins, 0., dTrange);
    TH2F channelsdt2("channelsdt2", "Segments dT vs SumEnergy", dTbins, 0., dTrange, hBins, 0., hRange);

    UShort_t Channel;
    UShort_t Energy;
    UShort_t Board;
    ULong64_t Timestamp;

    TNtuple *ntuple = (TNtuple *)f->Get("Data");
    ntuple->SetBranchAddress("Channel", &Channel);
    ntuple->SetBranchAddress("Energy", &Energy);
    ntuple->SetBranchAddress("Timestamp", &Timestamp);
    ntuple->SetBranchAddress("Board", &Board);

    Long64_t nentries = ntuple->GetEntries();

    ULong64_t tFired, eFired, dT;
    UShort_t chFired = 1000;

    for (Long64_t i = 0; i < nentries; i++) {
    	ntuple->GetEntry(i);
	if (Board <= 700){
    		if (Channel == 0 || Channel == 4) {
      			if (chFired > 10 || chFired == Channel) {
        			tFired = Timestamp;
        			eFired = Energy;
        			chFired = Channel;
      			} else {
        			dT = Timestamp - tFired;
        			chFired = 1000;
        			tFired = Timestamp;
        			channelsdt.Fill(dT);
        			channelsdt2.Fill(dT, eFired + Energy);
      			}
    		}
	}
    }

    f->cd();
    channelsdt.Write();
    channelsdt2.Write();

    f->Close();

    std::cout << "dT histograms written to " << filename << "_processed.root" << std::endl;
}

//--------------------------------------------------------------------------------------------//

void SumSpectra(TString filename = "./co_test_3_9_2024_0cm/UNFILTERED/SData_co_test_3_9_2024_0cm") {
  TFile *f = new TFile(filename + "_processed.root", "UPDATE");

  Float_t hRange = 15000;
  Long64_t dt = 20000.;
  Int_t hBins = 1024;

  Float_t Channel;
  Float_t CalibratedEnergy;
  Float_t Board;
  Float_t Timestamp;

  TNtuple *ntuple = (TNtuple *)f->Get("CalibratedData");
  ntuple->SetBranchAddress("Channel", &Channel);
  ntuple->SetBranchAddress("CalibratedEnergy", &CalibratedEnergy);
  ntuple->SetBranchAddress("Timestamp", &Timestamp);
  ntuple->SetBranchAddress("Board", &Board);

  Long64_t nentries = ntuple->GetEntries();
 
  TH2F sumH2("Sum", "Sum of segments 2D", hBins / 2., 0., hRange, 10, 0, 10);
  TH1I sumH("sumH", "Sum of segments", hBins, 0., hRange);

  Float_t tFired = 0, eFired = 0;
  Float_t chFired = 1000;
  bool save = false;
  Float_t mult = 0;

  for (Long64_t i = 0; i < nentries; i++) {
    ntuple->GetEntry(i);
    
    if (Board < 700) {  // Skip virtual channels
       if (((Timestamp - tFired) > dt && chFired != 1000 && tFired != 0) ||
        chFired == Channel) {
      if (chFired == Channel)
        mult = 1;
      else
        mult = 2;
      sumH2.Fill(eFired, mult);
      sumH.Fill(eFired);
      save = true;
    }

    if (save) {
      tFired = 0;
      eFired = 0;
      chFired = 1000;
      save = false;
      mult = 0;
    }

    if ((Channel == 0 || Channel == 4) && chFired == 1000) {
      tFired = Timestamp;
      eFired = CalibratedEnergy;
      chFired = Channel;
    } else if (Channel == 0 || Channel == 4)
      eFired += CalibratedEnergy;
  }
}
  sumH2.Write();
  sumH.Write();
  
  std::cout << "Sum spectra written to: " << filename << "_processed.root" << std::endl;

  f->Close();
}

//--------------------------------------------------------------------------------------------//

void MakeTVFiles(TString filename = "./co_test_3_9_2024_0cm/UNFILTERED/SData_co_test_3_9_2024_0cm") {
  TFile *f = new TFile(filename + "_processed.root", "READ");
  TH1I *sumAllh = (TH1I *)f->Get("sumH");
  TH1I *sumUp = (TH1I *)f->Get("chUp");
  TH1I *sumDown = (TH1I *)f->Get("chDown");
  TH1I *virtualCh = (TH1I *)f->Get("chVirtual");
  
  Int_t hBins = 4096;

  FILE *tv_upper_segment = fopen(filename + "_upper_segment.tv", "w");
  FILE *tv_summed = fopen(filename + "_summed.tv", "w");
  FILE *tv_lower_segment = fopen(filename + "_lower_segment.tv", "w");
  FILE *tv_virtual_channel = fopen(filename + "_virtual_channel.tv", "w");

  for (Long64_t i = 0; i < hBins - 1; i++) {
    fprintf(tv_upper_segment, "%10.3e\n", sumUp->GetBinContent(i));
  }

  for (Long64_t i = 0; i < hBins - 1; i++) {
    fprintf(tv_summed, "%10.3e\n", sumAllh->GetBinContent(i));
  }

  for (Long64_t i = 0; i < hBins - 1; i++) {
    fprintf(tv_lower_segment, "%10.3e\n", sumDown->GetBinContent(i));
  }
  
  for (Long64_t i = 0; i < hBins - 1; i++) {
    fprintf(tv_virtual_channel, "%10.3e\n", virtualCh->GetBinContent(i));
  }
}

TH1F* CalibrateHistogramX(TH1* hOriginal, Double_t gain = 1.0, Double_t offset = 0.0, TString newname = "hCalibrated") {
  Int_t nbins = hOriginal->GetNbinsX();
  Double_t xmin = hOriginal->GetXaxis()->GetXmin();
  Double_t xmax = hOriginal->GetXaxis()->GetXmax();

  // New calibrated range
  Double_t new_xmin = gain * xmin + offset;
  Double_t new_xmax = gain * xmax + offset;

  TH1F* hNew = new TH1F(newname, hOriginal->GetTitle(), nbins, new_xmin, new_xmax);
  //hNew->GetXaxis()->SetTitle("Calibrated Energy");

  for (Int_t i = 1; i <= nbins; ++i) {
    Double_t x = hOriginal->GetXaxis()->GetBinCenter(i);
    Double_t y = hOriginal->GetBinContent(i);
    Double_t x_calibrated = gain * x + offset;
    hNew->Fill(x_calibrated, y);
  }

  return hNew;
}

void PlotSpectraComparison(TString filename = "./co_test_3_9_2024_0cm/UNFILTERED/SData_co_test_3_9_2024_0cm") {
  TFile *f = new TFile(filename + "_processed.root", "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Error opening file for plotting." << std::endl;
    return;
  }

  TH1I *sumH = (TH1I *)f->Get("sumH");
  TH1I *chUp = (TH1I *)f->Get("chUp");
  TH1I *chDown = (TH1I *)f->Get("chDown");
  TH1I *chVirtual = (TH1I *)f->Get("chVirtual");

  if (!sumH || !chUp || !chDown || !chVirtual) {
    std::cerr << "One or more histograms missing in file." << std::endl;
    return;
  }

  // --- Apply X-axis calibration (adjust these values as needed)
  auto *hSum = CalibrateHistogramX(sumH, 1.0, 0.0, "hSum_cal");
  auto *hUp = CalibrateHistogramX(chUp, 1.0, 0.0, "hUp_cal");
  auto *hDown = CalibrateHistogramX(chDown, 1.0, 0.0, "hDown_cal");
  auto *hVirtual = CalibrateHistogramX(chVirtual, 1., 240.0, "hVirtual_cal");

  // Set line colors and styles
  hSum->SetLineColor(kBlack); hSum->SetLineWidth(2);
  hUp->SetLineColor(kBlue); //hUp->SetLineStyle(2);
  hDown->SetLineColor(kRed); //hDown->SetLineStyle(3);
  hVirtual->SetLineColor(kGreen + 2); //hVirtual->SetLineStyle(4);

  // Draw
  TCanvas *c = new TCanvas("cSpectra", "Spectra Comparison", 1000, 600);
  c->SetLogy();  // optional
  
  gStyle->SetOptStat(0);


  hSum->Draw("HIST");
  hUp->Draw("HIST SAME");
  hDown->Draw("HIST SAME");
  hVirtual->Draw("HIST SAME");

  // Add legend
  TLegend *leg = new TLegend(0.7, 0.65, 0.88, 0.88);
  leg->AddEntry(hSum, "Summed Segments", "l");
  leg->AddEntry(hUp, "Up Segment", "l");
  leg->AddEntry(hDown, "Down Segment", "l");
  leg->AddEntry(hVirtual, "Virtual Channel", "l");
  leg->Draw();

  c->Update();
  
  c->SaveAs(filename + "_comparison.png");
  std::cout << "Comparison plot saved as: " << filename << "_comparison.png" << std::endl;
}


//--------------------------------------------------------------------------------------------//

void Neoptolemos_sum_test(TString infile = "./co_test_3_9_2024_0cm/UNFILTERED/SData_co_test_3_9_2024_0cm") {
  ReadDataFile(infile);
  //CreateDTHistograms(infile);
  SumSpectra(infile);
  //MakeTVFiles(infile);
  //PlotSpectraComparison(infile);
  
}
