#include "Riostream.h"
#include "TBrowser.h"
#include "TFile.h"
#include "TH2.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TTree.h"
#include <stdio.h>
#include <string.h>

TString filename = "SDataR_SeperatePMTs";

Int_t hBins = 4096;

Long64_t dt = 500000;
// Long64_t dt = 1500000;

Int_t NbCh = 0;

TString filename_ext = ".root";

Float_t hRange = hBins - 1;

Float_t CoPeaks[2] = {1173., 1332.};
Double_t CalCoeffs[2] = {1, 0};
UShort_t maxCh = 0;

void ReadDataFile() {
  TFile *f = new TFile(filename + ".root", "READ");
  TFile *fnew = new TFile("NaI.root", "RECREATE");

  UShort_t Channel, Board, Energy;
  ULong64_t Timestamp;

  TNtuple *ntuple = (TNtuple *)f->Get("Data_R");
  ntuple->SetBranchAddress("Channel", &Channel);
  ntuple->SetBranchAddress("Board", &Board);
  ntuple->SetBranchAddress("Timestamp", &Timestamp);
  ntuple->SetBranchAddress("Energy", &Energy);

  TH1I ch0("ch0", "channel 0 spectrum", hBins, 0., hRange);
  TH1I ch1("ch1", "channel 1 spectrum", hBins, 0., hRange);
  TH1I ch2("ch2", "channel 2 spectrum", hBins, 0., hRange);
  TH1I ch3("ch3", "channel 3 spectrum", hBins, 0., hRange);

  Long64_t nentries = ntuple->GetEntries();
  Double_t mult_coeff, add_coeff;
  for (Long64_t i = 0; i < nentries; i++) {
    ntuple->GetEntry(i);
    maxCh = std::max(maxCh, Energy);
    if (Channel == 0)
      ch0.Fill(Energy);
    else if (Channel == 1)
      ch1.Fill(Energy);
    else if (Channel == 2)
      ch2.Fill(Energy);
    else if (Channel == 3)
      ch3.Fill(Energy);
  }

  NbCh = 4;

  f->Close();

  fnew->Write();
  fnew->Close();
}

void GainMatch() {
  TFile *f = new TFile(filename + ".root", "READ");
  TFile *fnew = new TFile("NaI.root", "UPDATE");
  Double_t chpeaks[2 * NbCh];
  Double_t maxchpeaks[2] = {0, 0};

  for (Int_t i = 0; i < NbCh; i++) {
    char nbstr[10];
    sprintf(nbstr, "ch%d", i);

    TH1F *hist = (TH1F *)fnew->Get(nbstr);
    cout << "GainMatching histogramm " << nbstr << " with name "
         << hist->GetName() << endl;
    Int_t dum_int, ch;
    Double_t ener, time;

    TSpectrum *s = new TSpectrum();
    Double_t *xpeaks = s->GetPositionX();
    //      Double_t *ypeaks = s->GetPositionY();
    s->Search(hist, .5, "nodraw");
    chpeaks[2 * i + 0] =
        (CoPeaks[0] + CoPeaks[1]) / (2 * (xpeaks[1] - xpeaks[0]));
    chpeaks[2 * i + 1] = ((CoPeaks[0] + CoPeaks[1]) / 2) *
                         (1 - (xpeaks[0] / (xpeaks[1] - xpeaks[0])));
    //      chpeaks[2 * i + 0] = xpeaks[0];
    //      chpeaks[2 * i + 1] = xpeaks[1];
    //      if (maxchpeaks[0] < xpeaks[0])
    //         for (size_t j = 0; j < 2; ++j)
    //            maxchpeaks[j] = xpeaks[j];
    cout << "For ch" << i << " I have found peaks at : " << chpeaks[2 * i + 0]
         << " and " << chpeaks[2 * i + 1] << endl;
  }

  //   CalCoeffs[0] = (CoPeaks[0] + CoPeaks[1]) / (2 * (maxchpeaks[1] -
  //   maxchpeaks[0])); CalCoeffs[1] = (CoPeaks[0] + CoPeaks[1]) - (CalCoeffs[0]
  //   * maxchpeaks[1]); hRange = 55000.;

  Int_t BoardNb, ChNb, Ener;
  Double_t GMEner;
  Long64_t TimeTAG;
  TNtuple *GMntuple = new TNtuple();
  GMntuple->SetName("NaIGM");
  GMntuple->SetTitle("data from ascii file Gain Matched");
  GMntuple->Branch("Channel", &ChNb, "ch/I");
  GMntuple->Branch("EnergyCh", &GMEner, "e/D");
  GMntuple->Branch("TimeTAG", &TimeTAG, "t/L");

  TH1I ch0c("ch0GM", "channel 0 spectrum Gain Matched", hBins, 0., hRange);
  TH1I ch1c("ch1GM", "channel 1 spectrum Gain Matched", hBins, 0., hRange);
  TH1I ch2c("ch2GM", "channel 2 spectrum Gain Matched", hBins, 0., hRange);
  TH1I ch3c("ch3GM", "channel 3 spectrum Gain Matched", hBins, 0., hRange);

  Double_t dTrange = 1e7;
  Long64_t dTbins = dTrange / 1000;
  TH1I ch01dt("ch01dt", "dT channel 0-1 spectrum", dTbins, 0., dTrange);
  TH1I ch23dt("ch23dt", "dT channel 2-3 spectrum", dTbins, 0., dTrange);
  TH2F ch01dt2("ch01dt2", "dT channel 0-1 spectrum", dTbins, 0., dTrange, hBins,
               0., hRange);
  TH2F ch23dt2("ch23dt2", "dT channel 2-3 spectrum", dTbins, 0., dTrange, hBins,
               0., hRange);

  UShort_t Channel, Board, Energy;
  ULong64_t Timestamp;
  TNtuple *ntuple = (TNtuple *)f->Get("Data_R");
  ntuple->SetBranchAddress("Channel", &Channel);
  ntuple->SetBranchAddress("Energy", &Energy);
  ntuple->SetBranchAddress("Timestamp", &Timestamp);

  Long64_t nentries = ntuple->GetEntries();
  Double_t mult_coeff, add_coeff;

  ULong64_t t01Fired, e01Fired, dT01, t23Fired, e23Fired, dT23;
  UShort_t ch01Fired = 1000;
  UShort_t ch23Fired = 1000;
  for (Long64_t i = 0; i < nentries; i++) {
    ntuple->GetEntry(i);
    //      mult_coeff = (maxchpeaks[0] - maxchpeaks[1]) / (chpeaks[2 * Channel
    //      + 0] - chpeaks[2 * Channel + 1]); add_coeff = (maxchpeaks[1] *
    //      chpeaks[2 * Channel + 0] - maxchpeaks[0] * chpeaks[2 * Channel + 1])
    //      / (chpeaks[2 * Channel + 0] - chpeaks[2 * Channel + 1]); GMEner =
    //      (Energy * mult_coeff) + add_coeff; GMEner = (GMEner * CalCoeffs[0])
    //      + CalCoeffs[1];
    GMEner = (Energy * chpeaks[2 * Channel + 0]) + chpeaks[2 * Channel + 1];
    GMEner = Energy;

    if (Channel == 0 || Channel == 1) {
      if (ch01Fired > 10 || ch01Fired == Channel) {
        t01Fired = Timestamp;
        e01Fired = GMEner;
        ch01Fired = Channel;
      } else {
        dT01 = Timestamp - t01Fired;
        ch01Fired = 1000;
        t01Fired = Timestamp;
        //            e01Fired = GMEner;
        //            ch01Fired = Channel;
        ch01dt.Fill(dT01);
        ch01dt2.Fill(dT01, e01Fired + GMEner);
      }
    } else {
      if (ch23Fired > 10 || ch23Fired == Channel) {
        t23Fired = Timestamp;
        e23Fired = GMEner;
        ch23Fired = Channel;
      } else {
        dT23 = Timestamp - t23Fired;
        ch23Fired = 1000;
        e23Fired = GMEner;
        //            t23Fired = Timestamp;
        //            ch23Fired = Channel;
        ch23dt.Fill(dT23);
        ch23dt2.Fill(dT23, e23Fired + GMEner);
      }
    }

    if (Channel == 0)
      ch0c.Fill(GMEner);
    else if (Channel == 1)
      ch1c.Fill(GMEner);
    else if (Channel == 2)
      ch2c.Fill(GMEner);
    else if (Channel == 3)
      ch3c.Fill(GMEner);
    ChNb = Channel;
    TimeTAG = Timestamp;
    GMntuple->Fill(ChNb, GMEner, TimeTAG);
  }
  ch0c.Write();
  ch1c.Write();
  ch2c.Write();
  ch3c.Write();

  ch01dt.Write();
  ch23dt.Write();
  ch01dt2.Write();
  ch23dt2.Write();

  GMntuple->Write();

  f->Close();
  fnew->Close();
}

void SumSpectra() {
  TFile *f = new TFile("NaI.root", "UPDATE");

  Int_t ChNb;
  Double_t Ener;
  Long64_t TimeTAG;

  //   TNtuple *ntuple = (TNtuple*)f->Get("NaI");
  TNtuple *ntuple = (TNtuple *)f->Get("NaIGM");
  ntuple->SetBranchAddress("Channel", &ChNb);
  ntuple->SetBranchAddress("EnergyCh", &Ener);
  ntuple->SetBranchAddress("TimeTAG", &TimeTAG);

  Long64_t nentries = ntuple->GetEntries();

  Int_t multSum01 = 0, multSum23 = 0, multSumAll = 0;
  Int_t eSum01 = 0, eSum23 = 0, eSumAll = 0;
  Long64_t tSum01, tSum23, tSumAll = -100;

  TH2F sum01h("Sum01", "Sum of channels 0 and 1", hBins, 0., hRange, 10, 0, 10);
  TH2F sum23h("Sum23", "Sum of channels 2 and 3", hBins, 0., hRange, 10, 0, 10);
  TH2F sumAllh("SumAll", "Sum of all channels", hBins, 0., hRange, 10, 0, 10);

  TH1F sum01h2("Sum01H", "Sum of channels 0 and 1 h", hBins, 0., hRange);
  TH1F sum23h2("Sum23H", "Sum of channels 2 and 3 h", hBins, 0., hRange);
  TH1F sumAllh2("SumAllH", "Sum of all channels h", hBins, 0., hRange);

  for (Long64_t i = 0; i < nentries; i++) {
    ntuple->GetEntry(i);
    Ener /= 2.;
    if (TimeTAG > tSumAll + dt) {
      if (multSum01 >= 1) {
        sum01h.Fill(eSum01, multSum01);
        sum01h2.Fill(eSum01);
      }
      if (multSum23 >= 1) {
        sum23h.Fill(eSum23, multSum23);
        sum23h2.Fill(eSum23);
      }
      if (multSumAll >= 1) {
        sumAllh.Fill(eSumAll, multSumAll);
        sumAllh2.Fill(eSumAll);
      }
      if (ChNb == 0 || ChNb == 1) {
        multSum01 = 1;
        multSum23 = 0;
        eSum01 = Ener;
        eSum23 = 0.;
        tSum01 = TimeTAG;
      } else if (ChNb == 2 || ChNb == 3) {
        multSum01 = 0;
        multSum23 = 1;
        eSum01 = 0.;
        eSum23 = Ener;
        tSum23 = TimeTAG;
      }
      multSumAll = 1;
      eSumAll = Ener;
      tSumAll = TimeTAG;
    } else {
      if (ChNb == 0 || ChNb == 1) {
        multSum01++;
        eSum01 += Ener;
      } else if (ChNb == 2 || ChNb == 3) {
        multSum23++;
        eSum23 += Ener;
      }
      multSumAll++;
      eSumAll += Ener;
    }
    if (i == nentries - 1) {
      if (multSum01 >= 1) {
        sum01h.Fill(eSum01, multSum01);
        sum01h2.Fill(eSum01);
      }
      if (multSum23 >= 1) {
        sum23h.Fill(eSum23, multSum23);
        sum23h2.Fill(eSum23);
      }
      if (multSumAll >= 1) {
        sumAllh.Fill(eSumAll, multSum01);
        sumAllh2.Fill(eSumAll);
      }
    }
  }

  sum01h.Write();
  sum23h.Write();
  sumAllh.Write();

  sum01h2.Write();
  sum23h2.Write();
  sumAllh2.Write();

  f->Close();
}

void NaI_sum() {
  ReadDataFile();
  GainMatch();
  SumSpectra();
}
