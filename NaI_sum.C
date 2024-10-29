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

Float_t CoPeaks[2] = {1173., 1332.};
Double_t CalCoeffs[2] = {1, 0};
UShort_t maxCh = 0;

void ReadDataFile(TString filename = "SDataR_SeperatePMTs", Int_t hBins = 4096,
                  TString filenameNew = "NaI") {
   TFile *f = new TFile(filename + ".root", "READ");
   TFile *fnew = new TFile(filenameNew + ".root", "RECREATE");

   Float_t hRange = hBins - 1;
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

void GainMatch(TString filename = "SDataR_SeperatePMTs", Int_t hBins = 4096,
               TString filenameNew = "NaI") {
   TFile *f = new TFile(filename + ".root", "READ");
   TFile *fnew = new TFile(filenameNew + ".root", "UPDATE");

   Float_t hRange = hBins - 1;
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
      Double_t sigma = 0.5;
      Double_t threashold = .05;
      s->Search(hist, sigma, "nodraw", threashold);
      while (s->GetNPeaks() > 2) {
         //         threashold *= 2;
         sigma *= 1.5;
         s->Search(hist, sigma, "nodraw", threashold);
         //         cout << "I have found: " << s->GetNPeaks() << " " << sigma << "
         //         " << threashold << endl;
      }
      Double_t *xpeaks = s->GetPositionX();
      //      for (Int_t kk = 0; kk < s->GetNPeaks(); ++kk)
      //         cout << "Peak " << kk << " " << xpeaks[kk] << endl;

      chpeaks[2 * i + 0] =
          (CoPeaks[0] + CoPeaks[1]) / (2 * (xpeaks[1] - xpeaks[0]));
      chpeaks[2 * i + 1] = (xpeaks[1] - (2 * xpeaks[0])) * chpeaks[2 * i + 0];

      cout << "For ch" << i << " I have found peaks at : " << chpeaks[2 * i + 0]
           << " and " << chpeaks[2 * i + 1] << endl;
   }

   hRange = 55000.;

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
      GMEner = (Energy * chpeaks[2 * Channel + 0]) + chpeaks[2 * Channel + 1];
      //      GMEner = Energy;

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

void SumSpectra(Int_t hBins = 4096, Long64_t dt = 100000,
                TString filenameNew = "NaI") {
   TFile *f = new TFile(filenameNew + ".root", "UPDATE");

   Float_t hRange = hBins - 1;
   hRange = 55000.;

   Int_t Channel;
   Double_t Ener;
   Long64_t Timestamp;

   //   TNtuple *ntuple = (TNtuple*)f->Get("NaI");
   TNtuple *ntuple = (TNtuple *)f->Get("NaIGM");
   ntuple->SetBranchAddress("Channel", &Channel);
   ntuple->SetBranchAddress("EnergyCh", &Ener);
   ntuple->SetBranchAddress("TimeTAG", &Timestamp);

   Long64_t nentries = ntuple->GetEntries();

   Int_t multSum01 = 0, multSum23 = 0, multSumAll = 0;
   Int_t eSum01 = 0, eSum23 = 0, eSumAll = 0;
   Long64_t tSum01, tSum23, tSumAll = -100;

   TH2F sum01h("Sum01", "Sum of channels 0 and 1", hBins, 0., hRange, 10, 0, 10);
   TH2F sum23h("Sum23", "Sum of channels 2 and 3", hBins, 0., hRange, 10, 0, 10);
   TH2F sumAllh("SumAll", "Sum of all channels", hBins, 0., hRange, 10, 0, 10);

   TH1F sum01h2("Sum01H", "Sum of channels 0 and 1 h", 2 * hBins, 0.,
                2 * hRange);
   TH1F sum23h2("Sum23H", "Sum of channels 2 and 3 h", 2 * hBins, 0.,
                2 * hRange);
   TH1F sumAllh2("SumAllH", "Sum of all channels h", 4 * hBins, 0., 4 * hRange);

   ULong64_t t01Fired = 0, e01Fired, dT01, t23Fired = 0, e23Fired, dT23;
   UShort_t ch01Fired = 1000;
   UShort_t ch23Fired = 1000;
   bool save01 = false;
   bool save23 = false;
   UShort_t mult01 = 0, mult23 = 0;

   for (Long64_t i = 0; i < nentries; i++) {
      ntuple->GetEntry(i);

      if (((Timestamp - t01Fired) > dt && ch01Fired != 1000 && t01Fired != 0) ||
          ch01Fired == Channel) {
         if (ch01Fired == Channel)
            mult01 = 1;
         else
            mult01 = 2;
         sum01h.Fill(e01Fired, mult01);
         sum01h2.Fill(e01Fired);
         save01 = true;
      }

      if (((Timestamp - t23Fired) > dt && ch23Fired != 1000 && t23Fired != 0) ||
          ch23Fired == Channel) {
         if (ch23Fired == Channel)
            mult23 = 1;
         else
            mult23 = 2;
         sum23h.Fill(e23Fired, mult23);
         sum23h2.Fill(e23Fired);
         save23 = true;
      }

      if (save01 || save23) {
         sumAllh2.Fill(e01Fired + e23Fired);
         sumAllh.Fill(e01Fired + e23Fired, mult01 + mult23);
         if (save01) {
            t01Fired = 0;
            e01Fired = 0;
            ch01Fired = 1000;
            save01 = false;
            mult01 = 0;
         }
         if (save23) {
            t23Fired = 0;
            e23Fired = 0;
            ch23Fired = 1000;
            save23 = false;
            mult23 = 0;
         }
      }

      if ((Channel == 0 || Channel == 1) && ch01Fired == 1000) {
         t01Fired = Timestamp;
         e01Fired = Ener;
         ch01Fired = Channel;
      } else if (Channel == 0 || Channel == 1)
         e01Fired += Ener;
      else if ((Channel == 2 || Channel == 3) && ch23Fired == 1000) {
         t23Fired = Timestamp;
         e23Fired = Ener;
         ch23Fired = Channel;
      } else if (Channel == 2 || Channel == 3)
         e23Fired += Ener;
   }

   sum01h.Write();
   sum23h.Write();
   sumAllh.Write();

   sum01h2.Write();
   sum23h2.Write();
   sumAllh2.Write();

   f->Close();
}

void MakeTVFiles(TString filename = "NaI", Int_t hBins = 4096,
                 TString filenameNew = "NaI") {
   TFile *f = new TFile(filenameNew + ".root", "READ");
   TH1F *sumAllh = (TH1F *)f->Get("SumAllH");
   TH1F *sum01h = (TH1F *)f->Get("Sum01H");
   TH1F *sum23h = (TH1F *)f->Get("Sum23H");

   FILE *tv_sum_all = fopen(filename + "_sum_all.tv", "w");
   FILE *tv_sum_01 = fopen(filename + "_sum_01.tv", "w");
   FILE *tv_sum_23 = fopen(filename + "_sum_23.tv", "w");

   for (Long64_t i = 0; i < hBins - 1; i++) {
      fprintf(tv_sum_all, "%10.3e\n", sumAllh->GetBinContent(i));
   }
   //   tv_sum_all->close();

   for (Long64_t i = 0; i < hBins - 1; i++) {
      fprintf(tv_sum_01, "%10.3e\n", sum01h->GetBinContent(i));
   }
   //   tv_sum_01->close();

   for (Long64_t i = 0; i < hBins - 1; i++) {
      fprintf(tv_sum_23, "%10.3e\n", sum23h->GetBinContent(i));
   }
   //   tv_sum_23->close();
   f->Close();

   filenameNew = "AllSumed/RAW/HcompassR_AllSumed_20230313_191539";
   f = new TFile(filenameNew + ".root", "READ");
   sumAllh = (TH1F *)f->Get("Energy/_R_EnergyCH0@DT5725SB_2131");
   tv_sum_all = fopen(filename + "_sum_all_orig.tv", "w");
   for (Long64_t i = 0; i < hBins - 1; i++) {
      fprintf(tv_sum_all, "%10.3e\n", sumAllh->GetBinContent(i));
   }
   //   tv_sum_all->close();
   f->Close();

   filenameNew = "SegmentSumed/RAW/HcompassR_SegmentSumed_20230313_191908";
   f = new TFile(filenameNew + ".root", "READ");
   sum01h = (TH1F *)f->Get("Energy/_R_EnergyCH0@DT5725SB_2131");
   sum23h = (TH1F *)f->Get("Energy/_R_EnergyCH1@DT5725SB_2131");
   tv_sum_01 = fopen(filename + "_sum_01_orig.tv", "w");
   tv_sum_23 = fopen(filename + "_sum_23_orig.tv", "w");
   for (Long64_t i = 0; i < hBins - 1; i++) {
      fprintf(tv_sum_01, "%10.3e\n", sum01h->GetBinContent(i));
      fprintf(tv_sum_23, "%10.3e\n", sum23h->GetBinContent(i));
   }
   //   tv_sum_01->close();
   //   tv_sum_23->close();
   f->Close();
}

void NaI_sum(TString infile = "SDataR_SeperatePMTs", Int_t bins = 4096,
             Long64_t delta_t = 15000, TString outfile = "NaI") {
   ReadDataFile(infile, bins, outfile);
   GainMatch(infile, bins, outfile);
   SumSpectra(bins, delta_t, outfile);
   MakeTVFiles(outfile, bins, outfile);
}
