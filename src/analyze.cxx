#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TLatex.h>
#include <TF1.h>

#include <iostream>
#include <sstream>
#include <string>

#include <math.h>

//#include "eps09_cxx/eps09.h"

#include <AtlasStyle.h>
#include <AtlasUtils.h>

#include <GlobalParams.h>
#include <Utilities.h>

using namespace atlashi;

string FormatCounts (int counts) {
  if (counts < 1000) return "";
  else if (1000 <= counts && counts < 10000) {
    string countsStr = FormatMeasurement (counts, 0, 1);
    countsStr = countsStr.substr(0, 1) + "k";
    return countsStr;
  }
  else if (10000 <= counts && counts < 100000) {
    string countsStr = FormatMeasurement (counts, 0, 2);
    countsStr = countsStr.substr(0, 2) + "k";
    return countsStr;
  }
  else if (100000 <= counts && counts < 1000000) {
    string countsStr = FormatMeasurement (counts, 0, 3);
    countsStr = countsStr.substr(0, 3) + "k";
    return countsStr;
  }
  else return "";
}


int GetNumInGroup1 (short iPtZ, short iCent) {
  if (iPtZ == 2) {
    if (iCent == 1)       return 260;
    else if (iCent == 2)  return 673;
    else if (iCent == 3)  return 782;
  }
  else if (iPtZ == 3) {
    if (iCent == 1)       return 148;
    else if (iCent == 2)  return 364;
    else if (iCent == 3)  return 443;
  }
  else if (iPtZ == 4) {
    if (iCent == 1)       return 52;
    else if (iCent == 2)  return 147;
    else if (iCent == 3)  return 145;
  }
  return 0;
}


int GetNumInGroup2 (short iPtZ, short iCent) {
  if (iPtZ == 2) {
    if (iCent == 1)       return 397;
    else if (iCent == 2)  return 828;
    else if (iCent == 3)  return 805;
  }
  else if (iPtZ == 3) {
    if (iCent == 1)       return 206;
    else if (iCent == 2)  return 484;
    else if (iCent == 3)  return 406;
  }
  else if (iPtZ == 4) {
    if (iCent == 1)       return 89;
    else if (iCent == 2)  return 171;
    else if (iCent == 3)  return 173;
  }
  return 0;
}


int main () {

  SetAtlasStyle();

  const int nSeeds = 401;
  TFile* inFile;
  TTree* inTree;

  //int code = 0;
  //int id1 = 0;
  //int id2 = 0;
  //float x1pdf = 0;
  //float x2pdf = 0;
  //float Q = 0;
  //bool isValence1 = false;
  //bool isValence2 = false;

  float z_pt = 0;
  float z_eta = 0;
  float z_phi = 0;
  float z_m = 0;
  int part_n = 0;
  float part_pt[10000];
  float part_eta[10000];
  float part_phi[10000];


  float trkPtYields[2][6];
  float trkPtYieldWeightsSq[2][6];
  float trkXYields[2][6];
  float trkXYieldWeightsSq[2][6];

  float trkPtYieldAvg[2][6];
  float trkPtYieldSumSq[2][6];
  float trkPtYieldCov[6];
  float trkXYieldAvg[2][6];
  float trkXYieldSumSq[2][6];
  float trkXYieldCov[6];

  TH1D* h_trk_pth_pull[3][3][6];
  TH1D* h_trk_xhz_pull[3][3][6];

  for (int iPtZ : {2, 3, 4}) {
    for (int iCent : {1, 2, 3}) {

      for (int iX = 0; iX < 6; iX++) {
        h_trk_pth_pull[iPtZ-2][iCent-1][iX] = new TH1D (Form ("h_trk_pth_pull_iPtZ%i_iCent%i_iPtTrk%i", iPtZ, iCent, iX), ";Pull;Counts", 12, -4, 4);
        h_trk_xhz_pull[iPtZ-2][iCent-1][iX] = new TH1D (Form ("h_trk_xhz_pull_iPtZ%i_iCent%i_iXHZ%i", iPtZ, iCent, iX), ";Pull;Counts", 12, -4, 4);
        h_trk_pth_pull[iPtZ-2][iCent-1][iX]->Sumw2 ();
        h_trk_xhz_pull[iPtZ-2][iCent-1][iX]->Sumw2 ();
      }

      const int nGroup1 = GetNumInGroup1 (iPtZ, iCent);
      const int nGroup2 = GetNumInGroup2 (iPtZ, iCent);

      for (int iSeed = 0; iSeed < nSeeds; iSeed++) {

        inFile = new TFile (Form ("output/seed%i_iPtZ%i_iCent%i.root", iSeed, iPtZ, iCent), "read");
        inTree = (TTree*) inFile->Get ("tree");

        //inTree->SetBranchAddress ("code",       &code);
        //inTree->SetBranchAddress ("id1",        &id1);
        //inTree->SetBranchAddress ("id2",        &id2);
        //inTree->SetBranchAddress ("x1pdf",      &x1pdf);
        //inTree->SetBranchAddress ("x2pdf",      &x2pdf);
        //inTree->SetBranchAddress ("Q",          &Q);
        //inTree->SetBranchAddress ("isValence1", &isValence1);
        //inTree->SetBranchAddress ("isValence2", &isValence2);
        inTree->SetBranchAddress ("z_pt",       &z_pt);
        inTree->SetBranchAddress ("z_eta",      &z_eta);
        inTree->SetBranchAddress ("z_phi",      &z_phi);
        inTree->SetBranchAddress ("z_m",        &z_m);
        inTree->SetBranchAddress ("part_n",     &part_n);
        inTree->SetBranchAddress ("part_pt",    &part_pt);
        inTree->SetBranchAddress ("part_eta",   &part_eta);
        inTree->SetBranchAddress ("part_phi",   &part_phi);

        const int nEvents = inTree->GetEntries ();
        assert (nEvents == nGroup1+nGroup2);

        for (int iEvent = 0; iEvent < nEvents; iEvent++) {
          inTree->GetEntry (iEvent);

          const int iGroup = (iEvent < nGroup1 ? 0 : 1); // 0 corresponds to group 1, 1 to group 2

          float tracks_1to2GeV = 0;
          float tracks_2to4GeV = 0;
          float tracks_4to8GeV = 0;
          float tracks_8to15GeV = 0;
          float tracks_15to30GeV = 0;
          float tracks_30to60GeV = 0;
          float tracks_1_60to1_30 = 0;
          float tracks_1_30to1_15 = 0;
          float tracks_1_15to1_8 = 0;
          float tracks_1_8to1_4 = 0;
          float tracks_1_4to1_2 = 0;
          float tracks_1_2to1 = 0;
      
          for (int iPart = 0; iPart < part_n; iPart++) {
            if (fabs (part_eta[iPart]) > 2.5)
              continue;
            if (DeltaPhi (part_phi[iPart], z_phi) < 3*pi/4)
              continue;

            const float trkpt = part_pt[iPart];
            const float xhz = trkpt / z_pt;

            if (trkpt < 60) {
              if (30 <= trkpt) tracks_30to60GeV += 1.;
              else if (15 <= trkpt) tracks_15to30GeV += 1.;
              else if (8 <= trkpt) tracks_8to15GeV += 1.;
              else if (4 <= trkpt) tracks_4to8GeV += 1.;
              else if (2 <= trkpt) tracks_2to4GeV += 1.;
              else if (1 <= trkpt) tracks_1to2GeV += 1.;
            }

            if (1./60. <= xhz) {
              if (xhz <= 1./30.) tracks_1_60to1_30 += 1.;
              else if (xhz <= 1./15.) tracks_1_30to1_15 += 1.;
              else if (xhz <= 1./8.) tracks_1_15to1_8 += 1.;
              else if (xhz <= 1./4.) tracks_1_8to1_4 += 1.;
              else if (xhz <= 1./2.) tracks_1_4to1_2 += 1.;
              else if (xhz <= 1.) tracks_1_2to1 += 1.;
            }
          } // end loop over iPart

          trkPtYields[iGroup][0] += tracks_1to2GeV;
          trkPtYields[iGroup][1] += tracks_2to4GeV;
          trkPtYields[iGroup][2] += tracks_4to8GeV;
          trkPtYields[iGroup][3] += tracks_8to15GeV;
          trkPtYields[iGroup][4] += tracks_15to30GeV;
          trkPtYields[iGroup][5] += tracks_30to60GeV;
          trkPtYieldWeightsSq[iGroup][0] += pow (tracks_1to2GeV, 2);
          trkPtYieldWeightsSq[iGroup][1] += pow (tracks_2to4GeV, 2);
          trkPtYieldWeightsSq[iGroup][2] += pow (tracks_4to8GeV, 2);
          trkPtYieldWeightsSq[iGroup][3] += pow (tracks_8to15GeV, 2);
          trkPtYieldWeightsSq[iGroup][4] += pow (tracks_15to30GeV, 2);
          trkPtYieldWeightsSq[iGroup][5] += pow (tracks_30to60GeV, 2);

          trkXYields[iGroup][0] += tracks_1_60to1_30;
          trkXYields[iGroup][1] += tracks_1_30to1_15;
          trkXYields[iGroup][2] += tracks_1_15to1_8;
          trkXYields[iGroup][3] += tracks_1_8to1_4;
          trkXYields[iGroup][4] += tracks_1_4to1_2;
          trkXYields[iGroup][5] += tracks_1_2to1;
          trkXYieldWeightsSq[iGroup][0] += pow (tracks_1_60to1_30, 2);
          trkXYieldWeightsSq[iGroup][1] += pow (tracks_1_30to1_15, 2);
          trkXYieldWeightsSq[iGroup][2] += pow (tracks_1_15to1_8, 2);
          trkXYieldWeightsSq[iGroup][3] += pow (tracks_1_8to1_4, 2);
          trkXYieldWeightsSq[iGroup][4] += pow (tracks_1_4to1_2, 2);
          trkXYieldWeightsSq[iGroup][5] += pow (tracks_1_2to1, 2);
        } // end loop over iEvent

        inFile->Close ();
        SaferDelete (inFile);

        // calculate pull on per-Z yields
        float yield1, yield2, sigma1, sigma2;
        for (int iX = 0; iX < 6; iX++) {
          yield1 = trkPtYields[0][iX] / nGroup1;
          //sigma1 = sqrt (trkPtYields[0][iX]) / nGroup1;
          sigma1 = sqrt (trkPtYieldWeightsSq[0][iX]) / nGroup1;
          yield2 = trkPtYields[1][iX] / nGroup2;
          //sigma2 = sqrt (trkPtYields[1][iX]) / nGroup2;
          sigma2 = sqrt (trkPtYieldWeightsSq[1][iX]) / nGroup2;

          trkPtYieldAvg[0][iX] += yield1;
          trkPtYieldAvg[1][iX] += yield2;
          trkPtYieldSumSq[0][iX] += yield1*yield1;
          trkPtYieldSumSq[1][iX] += yield2*yield2;
          trkPtYieldCov[iX] += yield1 * yield2;

          const float pullPt = (yield1 - yield2) / sqrt (sigma1*sigma1 + sigma2*sigma2);
          h_trk_pth_pull[iPtZ-2][iCent-1][iX]->Fill (pullPt);

          yield1 = trkXYields[0][iX] / nGroup1;
          //sigma1 = sqrt (trkXYields[0][iX]) / nGroup1;
          sigma1 = sqrt (trkXYieldWeightsSq[0][iX]) / nGroup1;
          yield2 = trkXYields[1][iX] / nGroup2;
          //sigma2 = sqrt (trkXYields[1][iX]) / nGroup2;
          sigma2 = sqrt (trkXYieldWeightsSq[1][iX]) / nGroup2;

          trkXYieldAvg[0][iX] += yield1;
          trkXYieldAvg[1][iX] += yield2;
          trkXYieldSumSq[0][iX] += yield1*yield1;
          trkXYieldSumSq[1][iX] += yield2*yield2;
          trkXYieldCov[iX] += yield1 * yield2;

          const float pullX = (yield1 - yield2) / sqrt (sigma1*sigma1 + sigma2*sigma2);
          h_trk_xhz_pull[iPtZ-2][iCent-1][iX]->Fill (pullX);
        }

        for (int iX = 0; iX < 6; iX++) {
          trkPtYields[0][iX]  = 0.;
          trkPtYields[1][iX]  = 0.;
          trkPtYieldWeightsSq[0][iX]  = 0.;
          trkPtYieldWeightsSq[1][iX]  = 0.;
          trkXYields[0][iX]   = 0.;
          trkXYields[1][iX]   = 0.;
          trkXYieldWeightsSq[0][iX]   = 0.;
          trkXYieldWeightsSq[1][iX]   = 0.;
        } // end loop over iX

      } // end loop over seeds

      cout << "iPtZ = " << iPtZ << ", iCent = " << iCent << endl;
      float sigma1, sigma2;
      for (int iX = 0; iX < 6; iX++) {
        trkPtYieldAvg[0][iX] = trkPtYieldAvg[0][iX] / nSeeds;
        trkPtYieldAvg[1][iX] = trkPtYieldAvg[1][iX] / nSeeds;
        sigma1 = sqrt (trkPtYieldSumSq[0][iX] / nSeeds - trkPtYieldAvg[0][iX]);
        sigma2 = sqrt (trkPtYieldSumSq[1][iX] / nSeeds - trkPtYieldAvg[1][iX]);
        trkPtYieldCov[iX] = (trkPtYieldCov[iX] / nSeeds) - trkPtYieldAvg[0][iX]*trkPtYieldAvg[1][iX];
        cout << "Pt correlation coeff.: " << trkPtYieldCov[iX]/(sigma1*sigma2) << endl;

        trkXYieldAvg[0][iX] = trkXYieldAvg[0][iX] / nSeeds;
        trkXYieldAvg[1][iX] = trkXYieldAvg[1][iX] / nSeeds;
        sigma1 = sqrt (trkXYieldSumSq[0][iX] / nSeeds - trkXYieldAvg[0][iX]);
        sigma2 = sqrt (trkXYieldSumSq[1][iX] / nSeeds - trkXYieldAvg[1][iX]);
        trkXYieldCov[iX] = (trkXYieldCov[iX] / nSeeds) - trkXYieldAvg[0][iX]*trkXYieldAvg[1][iX];

        //cout << "X covariance: " << trkPtYieldCov[iX] << endl;
      } // end loop over iX

      for (int iX = 0; iX < 6; iX++) {
        trkPtYieldAvg[0][iX] = 0.;
        trkPtYieldAvg[1][iX] = 0.;
        trkPtYieldSumSq[0][iX] = 0.;
        trkPtYieldSumSq[1][iX] = 0.;
        trkPtYieldCov[iX] = 0.;
        trkXYieldAvg[0][iX] = 0.;
        trkXYieldAvg[1][iX] = 0.;
        trkXYieldSumSq[0][iX] = 0.;
        trkXYieldSumSq[1][iX] = 0.;
        trkXYieldCov[iX] = 0.;
      } // end loop over iX
        

    } // end loop over iCent
  } // end loop over iPtZ


  TCanvas* c_pull_pth_iPtZ2 = new TCanvas ("c_pull_pth_iPtZ2", "", 800, 600);
  c_pull_pth_iPtZ2->Divide (4, 3);
  TCanvas* c_pull_xhz_iPtZ2 = new TCanvas ("c_pull_xhz_iPtZ2", "", 800, 600);
  c_pull_xhz_iPtZ2->Divide (4, 3);
  TCanvas* c_pull_pth_iPtZ3 = new TCanvas ("c_pull_pth_iPtZ3", "", 1000, 600);
  c_pull_pth_iPtZ3->Divide (5, 3);
  TCanvas* c_pull_xhz_iPtZ3 = new TCanvas ("c_pull_xhz_iPtZ3", "", 1000, 600);
  c_pull_xhz_iPtZ3->Divide (5, 3);
  TCanvas* c_pull_pth_iPtZ4 = new TCanvas ("c_pull_pth_iPtZ4", "", 1200, 600);
  c_pull_pth_iPtZ4->Divide (6, 3);
  TCanvas* c_pull_xhz_iPtZ4 = new TCanvas ("c_pull_xhz_iPtZ4", "", 1200, 600);
  c_pull_xhz_iPtZ4->Divide (6, 3);

  for (int iX = 0; iX < 4; iX++) {
    for (int iCent : {1, 2, 3}) {
      c_pull_pth_iPtZ2->cd (4*(iCent-1)+iX+1);
      h_trk_pth_pull[0][iCent-1][iX]->Draw ("hist");

      TF1* f = new TF1 (Form ("f_trk_pth_pull_iPtZ2_iCent%i_iX%i", iCent, iX), "gaus(0)", -4, 4);
      h_trk_pth_pull[0][iCent-1][iX]->Fit (f, "RN0Q");
      f->Draw ("same");

      float chi2 = f->GetChisquare ();
      int ndf = f->GetNDF ();
      myText (0.2, 0.86, kBlack, Form ("#sigma = %.2f", f->GetParameter (2)), 0.06);
      myText (0.2, 0.8, kBlack, Form ("#chi^{2}/dof = %.2f/%i", chi2, ndf), 0.06);

      c_pull_xhz_iPtZ2->cd (4*(iCent-1)+iX+1);
      h_trk_xhz_pull[0][iCent-1][iX+2]->Draw ("hist");

      f = new TF1 (Form ("f_trk_xhz_pull_iPtZ2_iCent%i_iX%i", iCent, iX+2), "gaus(0)", -4, 4);
      h_trk_xhz_pull[0][iCent-1][iX+2]->Fit (f, "RN0Q");
      f->Draw ("same");

      chi2 = f->GetChisquare ();
      ndf = f->GetNDF ();
      myText (0.2, 0.86, kBlack, Form ("#sigma = %.2f", f->GetParameter (2)), 0.06);
      myText (0.2, 0.8, kBlack, Form ("#chi^{2}/dof = %.2f/%i", chi2, ndf), 0.06);
    } // end loop over iCent
  } // end loop over x-values
  for (int iX = 0; iX < 5; iX++) {
    for (int iCent : {1, 2, 3}) {
      c_pull_pth_iPtZ3->cd (5*(iCent-1)+iX+1);
      h_trk_pth_pull[1][iCent-1][iX]->Draw ("hist");

      TF1* f = new TF1 (Form ("f_trk_pth_pull_iPtZ3_iCent%i_iX%i", iCent, iX), "gaus(0)", -4, 4);
      h_trk_pth_pull[1][iCent-1][iX]->Fit (f, "RN0Q");
      f->Draw ("same");

      float chi2 = f->GetChisquare ();
      int ndf = f->GetNDF ();
      myText (0.2, 0.86, kBlack, Form ("#sigma = %.2f", f->GetParameter (2)), 0.06);
      myText (0.2, 0.8, kBlack, Form ("#chi^{2}/dof = %.2f/%i", chi2, ndf), 0.06);

      c_pull_xhz_iPtZ3->cd (5*(iCent-1)+iX+1);
      h_trk_xhz_pull[1][iCent-1][iX+1]->Draw ("hist");

      f = new TF1 (Form ("f_trk_xhz_pull_iPtZ3_iCent%i_iX%i", iCent, iX+1), "gaus(0)", -4, 4);
      h_trk_xhz_pull[1][iCent-1][iX+1]->Fit (f, "RN0Q");
      f->Draw ("same");

      chi2 = f->GetChisquare ();
      ndf = f->GetNDF ();
      myText (0.2, 0.86, kBlack, Form ("#sigma = %.2f", f->GetParameter (2)), 0.06);
      myText (0.2, 0.8, kBlack, Form ("#chi^{2}/dof = %.2f/%i", chi2, ndf), 0.06);
    } // end loop over iCent
  } // end loop over x-values
  for (int iX = 0; iX < 6; iX++) {
    for (int iCent : {1, 2, 3}) {
      c_pull_pth_iPtZ4->cd (6*(iCent-1)+iX+1);
      h_trk_pth_pull[2][iCent-1][iX]->Draw ("hist");

      TF1* f = new TF1 (Form ("f_trk_pth_pull_iPtZ4_iCent%i_iX%i", iCent, iX), "gaus(0)", -4, 4);
      h_trk_pth_pull[2][iCent-1][iX]->Fit (f, "RN0Q");
      f->Draw ("same");

      float chi2 = f->GetChisquare ();
      int ndf = f->GetNDF ();
      myText (0.2, 0.86, kBlack, Form ("#sigma = %.2f", f->GetParameter (2)), 0.06);
      myText (0.2, 0.8, kBlack, Form ("#chi^{2}/dof = %.2f/%i", chi2, ndf), 0.06);

      c_pull_xhz_iPtZ4->cd (6*(iCent-1)+iX+1);
      h_trk_xhz_pull[2][iCent-1][iX]->Draw ("hist");

      f = new TF1 (Form ("f_trk_xhz_pull_iPtZ4_iCent%i_iX%i", iCent, iX), "gaus(0)", -4, 4);
      h_trk_xhz_pull[2][iCent-1][iX]->Fit (f, "RN0Q");
      f->Draw ("same");

      chi2 = f->GetChisquare ();
      ndf = f->GetNDF ();
      myText (0.2, 0.86, kBlack, Form ("#sigma = %.2f", f->GetParameter (2)), 0.06);
      myText (0.2, 0.8, kBlack, Form ("#chi^{2}/dof = %.2f/%i", chi2, ndf), 0.06);
    } // end loop over iCent
  } // end loop over x-values

  c_pull_pth_iPtZ2->SaveAs ("Plots/pull_pth_iPtZ2.pdf");
  c_pull_xhz_iPtZ2->SaveAs ("Plots/pull_xhz_iPtZ2.pdf");
  c_pull_pth_iPtZ3->SaveAs ("Plots/pull_pth_iPtZ3.pdf");
  c_pull_xhz_iPtZ3->SaveAs ("Plots/pull_xhz_iPtZ3.pdf");
  c_pull_pth_iPtZ4->SaveAs ("Plots/pull_pth_iPtZ4.pdf");
  c_pull_xhz_iPtZ4->SaveAs ("Plots/pull_xhz_iPtZ4.pdf");



  return 0;
}

void analyze () {
  main ();
}
