#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <sstream>

#include "Pythia8/Pythia.h"

#include <Utilities.h>
#include <GlobalParams.h>

using namespace Pythia8;
using namespace atlashi;

Pythia pythia;

int main (int argc, char** argv) {

  if (argc != 7) {
    std::cout << " usage: ./zgen SEED PTHATMIN MINZPT MAXZPT NEVT FILENAMEOUT" << std::endl;
    return 0;
  }

  // get arguments
  const int seed = atoi (argv[1]);
  const float ptHatMin = atof (argv[2]);
  const float minZPt = atof (argv[3]);
  const float maxZPt = atof (argv[4]);
  //const float minTrkPt = atof (argv[5]);
  //const float maxTrkPt = atof (argv[6]);
  const int nEvents = atoi (argv[5]);
  const string outFileName = string (argv[6]);

  // Process selection. LHC initialization. Histogram.
  pythia.readString ("Beams:eCM = 5020.");

  pythia.readString("23:onMode = off");
  pythia.readString("23:onIfAny = 11 13");

  pythia.readString ("WeakZ0:gmZmode = 2"); // set to Z's
  pythia.readString ("WeakSingleBoson:ffbar2gmZ = off");       // code 221  f + fbar  -> Z
  pythia.readString ("WeakDoubleBoson:ffbar2gmZgmZ = off");    // code 231  f + fbar  -> Z + Z
  pythia.readString ("WeakDoubleBoson:ffbar2ZW = off");        // code 232  f + fbar  -> Z + W
  pythia.readString ("WeakBosonAndParton:qqbar2gmZg = on");    // code 241  q + qbar  -> Z + g
  pythia.readString ("WeakBosonAndParton:qg2gmZq = on");       // code 242  q + g     -> Z + q
  pythia.readString ("WeakBosonAndParton:ffbar2gmZgm = off");  // code 243  f + fbar  -> Z + gamma
  pythia.readString ("WeakBosonAndParton:fgm2gmZf = off");     // code 244  f + gamma -> Z + f

  pythia.readString (Form ("PhaseSpace:pTHatMin = %g", ptHatMin));
  pythia.readString ("PhaseSpace:mHatMin = 60");
  pythia.readString ("Random:setSeed = on");
  pythia.readString (Form ("Random:seed = %i", seed));

  pythia.init ();

  TFile* outFile = new TFile (outFileName.c_str (), "RECREATE");

  int b_code = 0;
  int b_id1 = 0;
  int b_id2 = 0;
  float b_x1pdf = 0;
  float b_x2pdf = 0;
  float b_Q = 0;
  bool b_isValence1 = false;
  bool b_isValence2 = false;

  float b_z_pt = 0;
  float b_z_eta = 0;
  float b_z_phi = 0;
  float b_z_m = 0;

  TTree* outTree = new TTree ("tree", "tree");

  outTree->Branch ("code",  &b_code,  "code/I");
  outTree->Branch ("id1",   &b_id1,   "id1/I");
  outTree->Branch ("id2",   &b_id2,   "id2/I");
  outTree->Branch ("x1pdf", &b_x1pdf, "x1pdf/F");
  outTree->Branch ("x2pdf", &b_x2pdf, "x2pdf/F");
  outTree->Branch ("Q",     &b_Q,     "Q/F");

  outTree->Branch ("z_pt",  &b_z_pt,  "z_pt/F");
  outTree->Branch ("z_eta", &b_z_eta, "z_eta/F");
  outTree->Branch ("z_phi", &b_z_phi, "z_phi/F");
  outTree->Branch ("z_m",   &b_z_m,   "z_m/F");

  TLorentzVector l1, l2, z;
  
  for (int iEvent = 0; iEvent < nEvents; iEvent++) {
    if (!pythia.next ())
      continue;

    bool foundZ = false;

    for (int i = 0; i < pythia.event.size (); i++) {

      if (!pythia.event[i].isFinal()) continue; // check if in final state

      if (abs (pythia.event[i].id ()) != 11 && abs (pythia.event[i].id ()) != 13) continue; // check if electron or muon, resp.

      if (pythia.event[i].pT () < 20)
        continue;
      if (fabs (pythia.event[i].eta ()) > 2.5)
        continue;

      l1.SetPtEtaPhiM (pythia.event[i].pT (), pythia.event[i].eta (), pythia.event[i].phi (), pythia.event[i].m ());

      for (int j = 0; j < i; j++) {

        if (!pythia.event[j].isFinal()) continue; // check if in final state

        if (pythia.event[i].id () != -pythia.event[j].id ()) continue; // check if anti-particle of first particle

        if (pythia.event[j].pT () < 20)
          continue;
        if (fabs (pythia.event[j].eta ()) > 2.5)
          continue;

        l2.SetPtEtaPhiM (pythia.event[j].pT (), pythia.event[j].eta (), pythia.event[j].phi (), pythia.event[j].m ());

        // reconstruct Z
        z = l1+l2;

        b_z_pt  = z.Pt ();
        if (minZPt != -1 && maxZPt != -1)
          if (b_z_pt < minZPt || maxZPt < b_z_pt) continue; // pT cut on Z bosons
        b_z_m   = z.M ();
        if (b_z_m < 60 || 150 < b_z_m) continue; // loose invariant mass cut to make sure these are from Z decays

        b_z_eta = z.Eta ();
        b_z_phi = z.Phi ();
        foundZ = true;
        break;
      }
    }

    if (!foundZ) {
      iEvent--;
      continue;
    }

    b_code = pythia.info.code ();
    b_id1 = pythia.info.id1pdf ();
    b_id2 = pythia.info.id2pdf ();
    b_x1pdf = pythia.info.x1pdf ();
    b_x2pdf = pythia.info.x2pdf ();
    b_Q =  pythia.info.QFac ();
    
    b_isValence1 = pythia.info.isValence1 ();
    b_isValence2 = pythia.info.isValence2 ();

    outTree->Fill();

    if (iEvent % (nEvents/100) == 0)
      std::cout << iEvent / (nEvents/100) << "\% done...\r" << std::flush;
  }

  pythia.stat();
  
  outFile->Write();
  outFile->Close();

  return 0;
}
