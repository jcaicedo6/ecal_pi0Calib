#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "TCanvas.h"
// for ifstream
#include <fstream>
#include <TRandom3.h>


using namespace std;

// run the script with "root -l pi0InvM.C"

void pi0InvM() {
    const Double_t z_calo = 8.14; // position of calorimeter from the target in m
    const Double_t z_target = 0.09;    // position of target
    //const Double_t z_origin = 0.0;
    //const Double_t vertex_z = z_target - z_origin;  // position of vertex, where the pi0 is created, right in the middle of the target

    TChain *ch = new TChain("T");

    string path = "/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign_lowcur_April28/rootfiles/";
    ifstream infile("lists_runs/runfiles_3173_to_3192.txt");
    string filename;

    while (getline(infile, filename)) {
        ch->Add((path + filename).c_str());
    }

    // Initialize the pointers before setting the branch address
    Double_t ecal_e[100]; // Adjust the size if necessary
    Double_t ecal_x[100];
    Double_t ecal_y[100];

    ch->SetBranchAddress("earm.ecal.clus.e", &ecal_e);
    ch->SetBranchAddress("earm.ecal.clus.x", &ecal_x);
    ch->SetBranchAddress("earm.ecal.clus.y", &ecal_y);

    Double_t nclus = 0;
    ch->SetBranchAddress("earm.ecal.nclus", &nclus);

    Double_t clus_a_time[1000];
    ch->SetBranchAddress("earm.ecal.clus.atimeblk", &clus_a_time);

    Double_t clus_blk_atime[1000];
    Double_t clus_nblk[1000];
    Double_t clus_eblk[1000]; 
    ch->SetBranchAddress("earm.ecal.clus_blk.atime", &clus_blk_atime);
    ch->SetBranchAddress("earm.ecal.clus.nblk", &clus_nblk);  
    ch->SetBranchAddress("earm.ecal.clus.eblk", &clus_eblk);

    TH1F *h_pi0_mass = new TH1F("h_pi0_mass", "Pi0 Invariant Mass;M_{#pi^{0}} [GeV];Events", 80, 0, 0.6);
    TH2F *h_photons_time = new TH2F("h_photons_time", "Time of photons;earm.ecal.clus.atimeblk[0];earm.ecal.clus.atimeblk[1]", 100, 0, 540, 100, 0, 540);
    TH2F *h_xy = new TH2F("h_xy", "xy;earm.ecal.clus.x;earm.ecal.clus.y", 25, -1.5, 1.5, 25, -0.65, 0.6);
    TH1F *h_opening_angle = new TH1F("h_opening_angle", "Opening angle between #gamma#gamma [deg];Angle [deg];Events", 90, 0, 30);
    Long64_t nEvents = ch->GetEntries();
    cout << "Number of events: " << nEvents << endl;

    gRandom->SetSeed(0);
    for (Long64_t i = 0; i < nEvents; i++) {
        ch->GetEntry(i);
        // Check that there are at least two clusters
        //h_xy->Fill(ecal_x[0], ecal_y[0]);
        //h_xy->Fill(ecal_x[1], ecal_y[1]);
        if (nclus != 2) continue;
        if ((ecal_e[0] + ecal_e[1]) < 0.5) continue;
        //if (ecal_e[0] < 1) continue;
        //if (ecal_e[1] < 0.2) continue;
        if (ecal_e[0] < 0.2 || ecal_e[1] < 0.2) continue;    // Minimum cluster energy cut Tune the 0.2 GeV threshold based on your noise level.
        if (clus_nblk[0] < 2 || clus_nblk[1] < 2) continue;  // Minimum number of blocks per cluster
        
        //if (clus_nblk[0] < 2 || clus_nblk[1] < 2) continue; // Require at least 2 blocks per cluster
        //if (clus_eblk[0] < 0.05 || clus_eblk[1] < 0.05) continue; // Require significant energy in the highest-energy block

        Double_t deltaR = sqrt(pow(ecal_x[0] - ecal_x[1], 2) + pow(ecal_y[0] - ecal_y[1], 2));
        if (deltaR < 0.07) continue; // reject overlapping clusters

        //cout << "Processing event " << i << "\r";
        if (i % 1000 == 0) cout << "Processing event " << i << " / " << nEvents << endl;


        TVector3 pos1(ecal_x[0], ecal_y[0], z_calo);  // in m
        TVector3 pos2(ecal_x[1], ecal_y[1], z_calo);  // in m

        // uniform smearing ±15 cm around (0,0,0.09)
        double vertex_x_smeared = gRandom->Uniform(-0.15, +0.15);
        double vertex_y_smeared = gRandom->Uniform(-0.15, +0.15);
        double vertex_z_smeared = z_target + gRandom->Uniform(-0.15, +0.15);
        TVector3 vertex(vertex_x_smeared, vertex_y_smeared, vertex_z_smeared);
        TVector3 dir1 = (pos1 - vertex).Unit();
        TVector3 dir2 = (pos2 - vertex).Unit();

        TLorentzVector ph1(dir1.X() * ecal_e[0], dir1.Y() * ecal_e[0], dir1.Z() * ecal_e[0], ecal_e[0]);
        TLorentzVector ph2(dir2.X() * ecal_e[1], dir2.Y() * ecal_e[1], dir2.Z() * ecal_e[1], ecal_e[1]);


        Double_t pi0_mass = (ph1 + ph2).M();

        Double_t opening_angle = dir1.Angle(dir2) * (180.0 / TMath::Pi());
        if (opening_angle < 3.5 || opening_angle > 8) continue;  // The lower cut (e.g., < 6°) removes nearly collinear photon pairs → likely merged.
                                                                // The upper cut (e.g., > 80°) removes highly unphysical, possibly misreconstructed pairs.


        if ((ph1 + ph2).M() <= 0.4 && (ph1 + ph2).M() >= 0)
            { // cuts to ensure some quality and avoid bugs (inf, -nan)
                //cuts in time
                if (clus_a_time[0] > 100 && clus_a_time[0] < 300 && 
                    clus_a_time[1] > 100 && clus_a_time[1] < 300 && 
                    fabs(clus_a_time[0] - clus_a_time[1]) < 10) { // difference in time < 4 ns
                    h_pi0_mass->Fill(pi0_mass);
                    h_photons_time->Fill(clus_a_time[0], clus_a_time[1]);
                    h_opening_angle->Fill(opening_angle);
                    h_xy->Fill(ecal_x[0], ecal_y[0]);
                    h_xy->Fill(ecal_x[1], ecal_y[1]);       
                }
                
            }
           

        
    }

    TCanvas *c_pi0InvM = new TCanvas("c_pi0InvM", "c_pi0InvM", 800, 600);
    h_pi0_mass->Draw();

    TCanvas *c_photons_time = new TCanvas("c_photons_time", "c_photons_time", 800, 600);
    h_photons_time->Draw("colz");

    TCanvas *c_xy = new TCanvas("c_xy", "c_xy", 800, 600);
    h_xy->SetStats(0);
    h_xy->Draw("colz");

    TCanvas *c_opening = new TCanvas("c_opening", "Opening Angle", 800, 600);
    h_opening_angle->Draw();

    TFile *f = new TFile("pi0_mass.root", "RECREATE");
    h_pi0_mass->Write();
    f->Close();
}