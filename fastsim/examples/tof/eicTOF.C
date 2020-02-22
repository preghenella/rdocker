#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes);
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

std::map<int, uint> pidmap = { {11, 0}, {13, 1} , {211, 2}, {321, 3}, {2212, 4} };
const Char_t *pname[5]     = {"el", "mu", "pi", "ka", "pr"};
const Char_t *plabel[5]    = {"e^{#pm}", "#mu^{#pm}", "#pi^{#pm}", "K^{#pm}", "p+#bar{p}"};
Double_t pmass[5]          = {5.110e-04, 0.10566, 0.13957, 0.49368, 0.93827};

double c = 29.9792458;  // [cm/ns]
double tofR = 100.;     // [cm]
double tofZ = 200.;     // [cm]
double tofSigma = 0.02; // [ns]

double beta(Track *track)
{
  /** get info **/
  double tof  = track->TOuter * 1.e9; // [ns]
  double L    = track->L * 0.1; // [cm]
  return L / tof / c;
}

void pid(Track *track, double *deltat, double *nsigma, double etof)
{
  /** get info **/
  double tof  = track->TOuter * 1.e9; // [ns]
  double L    = track->L * 0.1;       // [cm]
  double p    = track->P;             // [GeV/c]
  double p2   = p * p;
  double ep   = p * 0.01;             // [GeV/c]
  double c    = 29.9792458;           // [cm/ns]
  double Lc   = L / c;

  /** perform PID **/
  for (Int_t ipart = 0; ipart < 5; ++ipart) {
    double mass2 = pmass[ipart] * pmass[ipart];
    double texp = Lc / p * TMath::Sqrt(mass2 + p2);
    double etexp = Lc * mass2 / p2 / TMath::Sqrt(mass2 + p2) * ep;    
    double sigma = TMath::Sqrt(etexp * etexp + etof * etof);
    deltat[ipart] = tof - texp;
    nsigma[ipart] = deltat[ipart] / sigma;
  }
}

void tof(std::string finname = "delphes.root",
	 std::string foutname = "tof.root")
{

  /** open file and connect to tree **/
  auto fin = TFile::Open(finname.c_str());
  auto tin = (TTree *)fin->Get("Delphes");
  auto nevents = tin->GetEntries();
  TClonesArray *events = new TClonesArray("HepMCEvent");
  tin->SetBranchAddress("Event", &events);
  TClonesArray *particles = new TClonesArray("GenParticle");
  tin->SetBranchAddress("particles", &particles);
  TClonesArray *tracks = new TClonesArray("Track");
  tin->SetBranchAddress("tracks", &tracks);

  /** TOF histograms **/
  auto hNhits = new TH1F("hNhits", "", 100, 0., 100.);
  auto hHitXY = new TH2F("hHitXY", ";x (cm);y (cm)", 500, -250., 250., 500, -250., 250.);
  auto hHitR  = new TH1F("hHitR", ";R (cm);", 250, 0., 250.);
  auto hHitZ  = new TH1F("hHitZ", ";z (cm);", 500, -250., 250.);
  auto hBetaP = new TH2F("hBetaP", ";#it{p} (GeV/#it{c});v/#it{c}", 1000, 0., 20., 1000, 0.1, 1.1);
  TH2 *hGenEtaP[5], *hAccEtaP[5], *hHitEtaP[5], *hDeltaT[5][5], *hNsigma[5][5];
  for (Int_t ipart = 0; ipart < 5; ++ipart) {
    hGenEtaP[ipart] = new TH2F(Form("hGenEtaP_%s", pname[ipart]), Form("All generated (%s);#eta; #it{p} (GeV/#it{c})", plabel[ipart]), 200, -5., 5., 200, 0., 10.);
    hAccEtaP[ipart] = new TH2F(Form("hAccEtaP_%s", pname[ipart]), Form("Generated in acceptance (%s);#eta; #it{p} (GeV/#it{c})", plabel[ipart]), 200, -5., 5., 200, 0., 10.);
    hHitEtaP[ipart] = new TH2F(Form("hHitEtaP_%s", pname[ipart]), Form("Hitting the detector (%s);#eta; #it{p} (GeV/#it{c})", plabel[ipart]), 200, -5., 5., 200, 0., 10.);
    
    for (Int_t iipart = 0; iipart < 5; ++iipart) {    
      hDeltaT[ipart][iipart] = new TH2F(Form("hDeltaT_%s_%s", pname[ipart], pname[iipart]),
					Form("True PID (%s) ;#it{p} (GeV/#it{c});#Deltat_{%s} (ns)", plabel[ipart], plabel[iipart]),
					200, 0., 20., 2000, -2.44, 2.44);
      hNsigma[ipart][iipart] = new TH2F(Form("hNsigma_%s_%s", pname[ipart], pname[iipart]),
					Form("True PID (%s);#it{p} (GeV/#it{c});n#sigma_{%s}", plabel[ipart], plabel[iipart]),
					200, 0., 20., 500, -25., 25.);
    }
  }

  /** PID variables **/
  Double_t deltat[5], nsigma[5];
  
  /** loop over events **/
  for (Int_t iev = 0; iev < nevents; ++iev) {

    tin->GetEntry(iev);

    /** loop over generated particles **/
    auto nparticles = particles->GetEntries();
    for (Int_t iparticle = 0; iparticle < nparticles; ++iparticle) {
      auto particle = (GenParticle *)particles->At(iparticle);
      if (!particle) continue;
      /** get information **/
      auto p_true          = particle->P;
      auto eta_true        = particle->Eta;
      auto pid_true        = particle->PID;
      auto is_primary_true = particle->T < 1.E-12;
      auto is_pikapr_true  = pidmap.count(abs(pid_true)) != 0;
      auto in_acceptance   = true;

      /** skip if not true primary pi/K/pr **/
      if (!is_primary_true || !is_pikapr_true) continue;
      auto ipart = pidmap[abs(pid_true)];

      /** fill eta-p map for all primary pi/K/pt generated **/
      hGenEtaP[ipart]->Fill(eta_true, p_true);
      
      /** skip if not in acceptance **/
      if (!in_acceptance) continue;
      
      /** fill eta-p map for primary pi/K/pt generated in acceptance **/
      hAccEtaP[ipart]->Fill(eta_true, p_true);
      
    } /** end of loop over generated particles **/
    
    /** loop over tracks **/
    auto ntracks = tracks->GetEntries();
    auto nhits = 0;
    for (Int_t itrk = 0; itrk < ntracks; ++itrk) {
      auto track = (Track *)tracks->At(itrk);
      if (!track) continue;
      /** get particle **/
      auto particle = (GenParticle *)track->Particle.GetObject();
      if (!particle) continue;
      /** get information **/
      auto p_true          = particle->P;
      auto eta_true        = particle->Eta;
      auto pid_true        = particle->PID;
      auto is_primary_true = particle->T < 1.E-12;
      auto is_pikapr_true  = pidmap.count(abs(pid_true)) != 0;
      auto x               = track->XOuter * 0.1; // [cm]
      auto y               = track->YOuter * 0.1; // [cm]
      auto z               = track->ZOuter * 0.1; // [cm]
      auto p               = track->P;            // [GeV]

      /** check track hits TOF **/
      auto hasTOF = fabs(hypot(x, y) - tofR) < 0.001 && fabs(z) < tofZ;
      if (!hasTOF) continue;
      
      nhits++;
      /** fill hit map with all tracks **/
      hHitXY->Fill(x, y);
      hHitR->Fill(hypot(x, y));
      hHitZ->Fill(z);

      /** skip if not true primary pi/K/pr **/
      if (!is_primary_true || !is_pikapr_true) continue;
      auto ipart = pidmap[abs(pid_true)];

      /** fill eta-p map for primary pi/K/pt actually arriving in acceptance **/
      hHitEtaP[ipart]->Fill(eta_true, p_true);

      /** check if has cherenkov light **/
      //      if (cherenkov(track)) continue;
      
      /** fill beta-p map for primary pi/K/pt actually arriving in acceptance **/
      hBetaP->Fill(p, beta(track));

      /** perform particle identification **/
      pid(track, deltat, nsigma, tofSigma);

      for (Int_t iipart = 0; iipart < 5; ++iipart) {    
	hDeltaT[ipart][iipart]->Fill(p, deltat[iipart]);
	hNsigma[ipart][iipart]->Fill(p, nsigma[iipart]);
      }
      
    } /** end of loop over TOF tracks **/

    hNhits->Fill(nhits);

  } /** end of loop over events **/

  /** save output and close **/
  auto fout = TFile::Open(foutname.c_str(), "RECREATE");
  hNhits->Write();
  hHitXY->Write();
  hHitR->Write();
  hHitZ->Write();
  hBetaP->Write();
  for (Int_t ipart = 0; ipart < 5; ++ipart) {
    hGenEtaP[ipart]->Write();
    hAccEtaP[ipart]->Write();
    hHitEtaP[ipart]->Write();
    for (Int_t iipart = 0; iipart < 5; ++iipart) {    
      hDeltaT[ipart][iipart]->Write();
      hNsigma[ipart][iipart]->Write();
    }
  }
  fout->Close();
  fin->Close();
}

