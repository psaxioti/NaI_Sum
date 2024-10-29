#include "Riostream.h"
#include "TBrowser.h"
#include "TFile.h"
#include "TH2.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TTree.h"
#include <stdio.h>
#include <string.h>

Int_t NbCh = 0;

UShort_t maxCh = 0;
Double_t up_mult_coeff = 7.34;
Double_t up_add_coeff = -1979.21;
Double_t down_mult_coeff = 7.10;
Double_t down_add_coeff = -1817.78;

Int_t hBins = 4096;
Float_t hRange = 20000.0;

Long64_t delta_t = 200000;

//--------------------------------------------------------------------------------------------//

void SetUp(Double_t up_mult = 7.34, Double_t up_add = -1979.21,
           Double_t down_mult = 7.10, Double_t down_add = -1817.78,
           Int_t bins = 4096, Float_t range = 20000, Long64_t dt = 200000) {
  up_mult_coeff = up_mult;
  up_add_coeff = up_add;
  down_mult_coeff = down_mult;
  down_add_coeff = down_add;

  hBins = bins;
  hRange = range;

  delta_t = dt;
}

//--------------------------------------------------------------------------------------------//

void ReadDataFile(TString filename =
                      "./240910_ThinAl_Scale_6A_998keV_p/FILTERED/"
                      "SDataF_240910_ThinAl_Scale_6A_998keV_p.root",
                  TString filenameNew = "NaI") {
  TFile *f = new TFile(filename + ".root", "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  } else {
    cout << "Reading " << filename << " file" << endl;
  }
  TFile *fnew = new TFile(filenameNew + ".root", "RECREATE");

  Float_t hRange = 4096;
  UShort_t Channel, Board, Energy;
  ULong64_t Timestamp;

  TNtuple *ntuple = (TNtuple *)f->Get("Data_F");
  ntuple->SetBranchAddress("Channel", &Channel);
  ntuple->SetBranchAddress("Board", &Board);
  ntuple->SetBranchAddress("Timestamp", &Timestamp);
  ntuple->SetBranchAddress("Energy", &Energy);

  if (!ntuple) {
    std::cerr << "TNtuple 'Data' not found in file." << std::endl;
    return;
  }

  TH1I chUp("chUp", "Half Up spectrum", hBins, 0., hRange);
  chUp.GetXaxis()->SetTitle("Channels");
  TH1I chDown("chDown", "Half Down spectrum", hBins, 0., hRange);
  chDown.GetXaxis()->SetTitle("Channels");

  Long64_t nentries = ntuple->GetEntries();
  for (Long64_t i = 0; i < nentries; i++) {
    ntuple->GetEntry(i);
    maxCh = std::max(maxCh, Energy);
    if (Channel == 1)
      chUp.Fill(Energy);
    else if (Channel == 4)
      chDown.Fill(Energy);
  }
  cout << "Original histograms created" << endl;
  NbCh = 2;

  f->Close();

  fnew->Write();
  fnew->Close();
}

//--------------------------------------------------------------------------------------------//

void Calibration(TString filename, TString filenameNew = "NaI") {
  TFile *f = new TFile(filename + ".root", "READ");
  TFile *fnew = new TFile(filenameNew + ".root", "UPDATE");

  TH1I chUpc("chUpCal", "Up Half spectrum Calibrated", hBins, 0., hRange);
  chUpc.GetXaxis()->SetTitle("Energy (keV)");
  TH1I chDownc("chDownCal", "Down Half spectrum Calibrated", hBins, 0., hRange);
  chDownc.GetXaxis()->SetTitle("Energy (keV)");

  UShort_t Channel, Energy;
  ULong64_t Timestamp;
  TNtuple *ntuple = (TNtuple *)f->Get("Data_F");
  ntuple->SetBranchAddress("Channel", &Channel);
  ntuple->SetBranchAddress("Energy", &Energy);
  ntuple->SetBranchAddress("Timestamp", &Timestamp);

  Long64_t nentries = ntuple->GetEntries();
  Int_t BoardNb, ChNb, Ener;
  Double_t CalEner;
  Long64_t TimeTAG;
  TNtuple *Calntuple = new TNtuple();
  Calntuple->SetName("NaICal");
  Calntuple->SetTitle("data Calibrated");
  Calntuple->Branch("Channel", &ChNb, "ch/I");
  Calntuple->Branch("EnergyCh", &CalEner, "e/D");
  Calntuple->Branch("TimeTAG", &TimeTAG, "t/L");

  for (Long64_t i = 0; i < nentries; i++) {
    ntuple->GetEntry(i);
    ChNb = Channel;
    TimeTAG = Timestamp;
    if (Channel == 1) {
      CalEner = Energy * up_mult_coeff + up_add_coeff;
      chUpc.Fill(CalEner);
    } else if (Channel == 4) {
      CalEner = Energy * down_mult_coeff + down_add_coeff;
      chDownc.Fill(CalEner);
    } else {
      CalEner = Energy;
    }

    if (Calntuple) {
      Calntuple->Fill(ChNb, CalEner, TimeTAG);
    } else {
      std::cerr << "TNtuple 'NaICal' not created properly." << std::endl;
    }
  }
  cout << "Calibrated histograms created" << endl;
  chUpc.Write();
  chDownc.Write();
  Calntuple->Write();

  f->Close();
  fnew->Close();
}

//--------------------------------------------------------------------------------------------//

void CreateDTHistograms(TString filename =
                            "./240910_ThinAl_Scale_6A_998keV_p/FILTERED/"
                            "SDataF_240910_ThinAl_Scale_6A_998keV_p",
                        TString filenameNew = "NaI") {

  TFile *f = new TFile(filename + ".root", "READ");
  TFile *fnew = new TFile(filenameNew + ".root", "UPDATE");

  Double_t dTrange = 1e6;
  Long64_t dTbins = dTrange / 1000;
  TH1I chdt("chdt", "Segments dT spectrum", dTbins, 0., dTrange);
  chdt.GetXaxis()->SetTitle("dT (ns)");

  TH2F chdt2("chdt2", "Segments dT spectrum", dTbins, 0., dTrange, hBins, 0.,
             hRange);

  UShort_t Channel, Board, Energy;
  ULong64_t Timestamp;
  TNtuple *ntuple = (TNtuple *)f->Get("Data");
  ntuple->SetBranchAddress("Channel", &Channel);
  ntuple->SetBranchAddress("Energy", &Energy);
  ntuple->SetBranchAddress("Timestamp", &Timestamp);

  Long64_t nentries = ntuple->GetEntries();

  ULong64_t t01Fired, e01Fired, dT01;
  UShort_t chUp1Fired = 1000;

  for (Long64_t i = 0; i < nentries; i++) {
    ntuple->GetEntry(i);

    if (Channel == 1 || Channel == 4) {
      if (chUp1Fired > 10 || chUp1Fired == Channel) {
        t01Fired = Timestamp;
        e01Fired = Energy;
        chUp1Fired = Channel;
      } else {
        dT01 = Timestamp - t01Fired;
        chUp1Fired = 1000;
        t01Fired = Timestamp;
        chdt.Fill(dT01);
        chdt2.Fill(dT01, e01Fired + Energy);
      }
    }
  }
  cout << "Histogram of time difference between Upper & Lower Segment created"
       << endl;
  chdt.Write();
  chdt2.Write();

  f->Close();
  fnew->Close();
}

//--------------------------------------------------------------------------------------------//

void SumSpectra(TString filenameNew = "NaI") {
  TFile *f = new TFile(filenameNew + ".root", "UPDATE");

  Int_t Channel;
  Double_t Ener;
  Long64_t Timestamp;

  TNtuple *ntuple = (TNtuple *)f->Get("NaICal");
  ntuple->SetBranchAddress("Channel", &Channel);
  ntuple->SetBranchAddress("EnergyCh", &Ener);
  ntuple->SetBranchAddress("TimeTAG", &Timestamp);

  Long64_t nentries = ntuple->GetEntries();

  TH1I sumH("sumH", "Sum of segments", hBins, 0., hRange);
  sumH.GetXaxis()->SetTitle("Energy (keV)");
  TH1I Original_Channel1("Original_Channel1", "Original_Channel1", hBins, 0.,
                         hRange);
  Original_Channel1.GetXaxis()->SetTitle("Energy (keV)");
  TH1I Original_Channel2("Original_Channel2", "Original_Channel2", hBins, 0.,
                         hRange);
  Original_Channel2.GetXaxis()->SetTitle("Energy (keV)");

  ULong64_t time_fired = 0;      // entry's timestamp
  ULong64_t energy_fired = 0;    // entry's energy
  UShort_t channel_fired = 1000; // 1000 for the 1st the first entry or for
                                 // every time we sum energies
  bool save_info = false;
  UShort_t channel_info = 0; // Channel number from the summed energy
  ULong64_t energy_info = 0; // energy value from the summed energy

  for (Long64_t i = 0; i < nentries; i++) {
    ntuple->GetEntry(i);

    if (((Timestamp - time_fired) > delta_t && channel_fired != 1000 &&
         time_fired != 0) ||
        channel_fired == Channel) {

      sumH.Fill(energy_fired);
      // to do
      if ((energy_fired > 12800) && (energy_fired < 13300)) {
        if (channel_info == 1) {
          Original_Channel1.Fill(energy_info);
        } else if (channel_info == 4) {
          Original_Channel2.Fill(energy_info);
        }
      }
      save_info = true;
    }

    if (save_info) {
      time_fired = 0;
      energy_fired = 0;
      channel_fired = 1000;
      save_info = false;
    }

    if ((Channel == 1 || Channel == 4) && channel_fired == 1000) {
      time_fired = Timestamp;
      energy_fired = Ener;
      channel_fired = Channel;
    } else if (Channel == 1 || Channel == 4)
      energy_fired += Ener;

    channel_info = Channel;
    energy_info = Ener;
  }

  cout << "Histogram of two segments summed created" << endl;
  sumH.Write();
  Original_Channel1.Write();
  Original_Channel2.Write();

  f->Close();
}
//--------------------------------------------------------------------------------------------//

void MakeTVFiles(TString filename =
                     "./240910_ThinAl_Scale_6A_998keV_p/FILTERED/"
                     "SDataF_240910_ThinAl_Scale_6A_998keV_p",
                 Int_t hBins = 4096, TString filenameNew = "NaI") {
  TFile *f = new TFile(filenameNew + ".root", "READ");
  TH1I *sumAllh = (TH1I *)f->Get("sumH");
  TH1I *sumUp = (TH1I *)f->Get("chUp");
  TH1I *sumDown = (TH1I *)f->Get("chDown");

  FILE *tv_upper_segment = fopen(filename + "_upper_segment.tv", "w");
  FILE *tv_summed = fopen(filename + "_summed.tv", "w");
  FILE *tv_lower_segment = fopen(filename + "_lower_segment.tv", "w");

  for (Long64_t i = 0; i < hBins - 1; i++) {
    fprintf(tv_upper_segment, "%10.3e\n", sumUp->GetBinContent(i));
  }

  for (Long64_t i = 0; i < hBins - 1; i++) {
    fprintf(tv_summed, "%10.3e\n", sumAllh->GetBinContent(i));
  }

  for (Long64_t i = 0; i < hBins - 1; i++) {
    fprintf(tv_lower_segment, "%10.3e\n", sumDown->GetBinContent(i));
  }
}

//--------------------------------------------------------------------------------------------//

void NaI_sum(
    TString infile = "./240910_ThinAl_Scale_6A_998keV_p/FILTERED/"
                     "SDataF_240910_ThinAl_Scale_6A_998keV_p",
    TString outfile = "NaI") {
  SetUp(7.34, -1979.21, 7.10, -1817.78, 2048, 20000, 200000);
  ReadDataFile(infile, outfile);
  Calibration(infile, outfile);
  SumSpectra(outfile);
  //  MakeTVFiles(infile, bins, outfile);
  // CreateDTHistograms(infile, outfile);
}
