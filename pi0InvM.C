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
    const Double_t z_calo = 6; // position of calorimeter from the target in m
    const Double_t z_target = 0.09;    // position of target
    //const Double_t z_origin = 0.0;
    //const Double_t vertex_z = z_target - z_origin;  // position of vertex, where the pi0 is created, right in the middle of the target

    TChain *ch = new TChain("T");

    //string path = "/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign_lowcur_April28/rootfiles/";
    //ifstream infile("lists_runs/runfiles_3173_to_3192.txt");

    string path = "/adaqfs/home/a-onl/sbs/Rootfiles/";
    ifstream infile("lists_runs/runfiles_3637.txt");
    string filename;

    while (getline(infile, filename)) {
        ch->Add((path + filename).c_str());
    }

    // Initialize the pointers before setting the branch address
    Double_t ecal_e[100]; // Adjust the size if necessary
    Double_t ecal_x[100];
    Double_t ecal_y[100];

    ch->SetBranchStatus("*", 0); // disable all Branches
    ch->SetBranchAddress("earm.ecal.clus.e", &ecal_e);
    ch->SetBranchAddress("earm.ecal.clus.x", &ecal_x);
    ch->SetBranchAddress("earm.ecal.clus.y", &ecal_y);

    Double_t nclus;
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
    TH2F *h_photons_time = new TH2F("h_photons_time", "Time of photons;earm.ecal.clus.atimeblk[icl];earm.ecal.clus.atimeblk[jcl]", 100, 0, 540, 100, 0, 540);
    TH2F *h_xy = new TH2F("h_xy", "xy;earm.ecal.clus.x;earm.ecal.clus.y", 25, -1.5, 1.5, 25, -0.65, 0.6);
    TH1F *h_opening_angle = new TH1F("h_opening_angle", "Opening angle between #gamma#gamma [deg];Angle [deg];Events", 90, 0, 30);
    Long64_t nEvents = ch->GetEntries();
    cout << "Number of events: " << nEvents << endl;

    gRandom->SetSeed(0);
    for (Long64_t i = 0; i < nEvents; i++) {
        ch->GetEntry(i);
        // Check that there are at least two clusters
        // Print progress
        if (i % 10000 == 0) {
            cout << "Processing event " << i << " / " << nEvents << "\r";
            cout.flush(); // Ensure the output is immediately displayed
        }
        //if (nclus < 2 || nclus > 3) continue;
        if (nclus < 2) continue;
        // Best pair of cluster per event
        double best_dt = 1e9;
        int best_icl = -1, best_jcl = -1;
        for (int icl = 0; icl < nclus; ++icl) {
            for (int jcl = icl + 1; jcl < nclus; ++jcl) {
                if ((ecal_e[icl] + ecal_e[jcl]) < 0.5) continue;
                if (ecal_e[icl] < 0.2 || ecal_e[jcl] < 0.2) continue;
                if (clus_nblk[icl] < 2 || clus_nblk[jcl] < 2) continue;
                Double_t deltaR = sqrt(pow(ecal_x[icl] - ecal_x[jcl], 2) + pow(ecal_y[icl] - ecal_y[jcl], 2));
                if (deltaR < 0.07) continue;
                if (clus_a_time[icl] < 100 || clus_a_time[icl] > 300) continue;
                if (clus_a_time[jcl] < 100 || clus_a_time[jcl] > 300) continue;
                double dt = fabs(clus_a_time[icl] - clus_a_time[jcl]);
                if (dt < best_dt && dt < 4) { // time window 
                    best_dt = dt;
                    best_icl = icl;
                    best_jcl = jcl;
                }
            }
        }
        // Find indices of the two clusters with the highest energies
        int imax1 = -1, imax2 = -1;                   // Will store the indices of the highest and second-highest energy clusters
        double e1 = -1, e2 = -1;                    // Will store the values of the highest and second-highest energies
        for (int icl = 0; icl < nclus; ++icl) {
            if (ecal_e[icl] > e1) {
                // If this cluster has the highest energy so far:
                e2 = e1; imax2 = imax1;               // The old highest becomes the second-highest, The old highest index becomes the second-highest index
                e1 = ecal_e[icl]; imax1 = icl;       // Update highest energy and highest index
            } else if (ecal_e[icl] > e2) {
                // If this cluster is not the highest, but is higher than the second-highest:
                e2 = ecal_e[icl]; imax2 = icl;       // Update second-highest energy and the second-highest index
            }
        }
        if (best_icl >= 0 && best_jcl >= 0) {
            if (!((best_icl == imax1 && best_jcl == imax2) || (best_icl == imax2 && best_jcl == imax1))) {
                continue; // skip if not the two highest-energy clusters
            }
            TVector3 pos1(ecal_x[best_icl], ecal_y[best_icl], z_calo);  // in m
            TVector3 pos2(ecal_x[best_jcl], ecal_y[best_jcl], z_calo);  // in m

            TVector3 vertex(0, 0, z_target);
            TVector3 dir1 = (pos1 - vertex).Unit();
            TVector3 dir2 = (pos2 - vertex).Unit();

            TLorentzVector ph1(dir1.X() * ecal_e[best_icl], dir1.Y() * ecal_e[best_icl], dir1.Z() * ecal_e[best_icl], ecal_e[best_icl]);
            TLorentzVector ph2(dir2.X() * ecal_e[best_jcl], dir2.Y() * ecal_e[best_jcl], dir2.Z() * ecal_e[best_jcl], ecal_e[best_jcl]);


            Double_t pi0_mass = (ph1 + ph2).M();

            Double_t opening_angle = dir1.Angle(dir2) * (180.0 / TMath::Pi());
            if (opening_angle < 3.5 || opening_angle > 8) continue;  // The lower cut (e.g., < 6°) removes nearly collinear photon pairs → likely merged.
                                                                    // The upper cut (e.g., > 80°) removes highly unphysical, possibly misreconstructed pairs.
            if (pi0_mass < 0 || pi0_mass > 0.4) continue;
            h_pi0_mass->Fill(pi0_mass);
            h_photons_time->Fill(clus_a_time[best_icl], clus_a_time[best_jcl]);
            h_opening_angle->Fill(opening_angle);
            h_xy->Fill(ecal_x[best_icl], ecal_y[best_icl]);
            h_xy->Fill(ecal_x[best_jcl], ecal_y[best_jcl]);
        
        /*// Loop over all unique pairs of clusters
        for (int icl = 0; icl < nclus; ++icl) {
            for (int jcl = icl + 1; jcl < nclus; ++jcl) {
                if ((ecal_e[icl] + ecal_e[jcl]) < 0.5) continue;
                //if (ecal_e[icl] < 1) continue;
                //if (ecal_e[jcl] < 0.2) continue;
                if (ecal_e[icl] < 0.2 || ecal_e[jcl] < 0.2) continue;    // Minimum cluster energy cut Tune the 0.2 GeV threshold based on your noise level.
                if (clus_nblk[icl] < 2 || clus_nblk[jcl] < 2) continue;  // Minimum number of blocks per cluster
                
                //if (clus_nblk[icl] < 2 || clus_nblk[jcl] < 2) continue; // Require at least 2 blocks per cluster
                //if (clus_eblk[icl] < 0.05 || clus_eblk[jcl] < 0.05) continue; // Require significant energy in the highest-energy block

                Double_t deltaR = sqrt(pow(ecal_x[icl] - ecal_x[jcl], 2) + pow(ecal_y[icl] - ecal_y[jcl], 2));
                if (deltaR < 0.07) continue; // reject overlapping clusters

                TVector3 pos1(ecal_x[icl], ecal_y[icl], z_calo);  // in m
                TVector3 pos2(ecal_x[jcl], ecal_y[jcl], z_calo);  // in m

                // uniform smearing ±15 cm around (0,0,0.09)
                //double vertex_x_smeared = gRandom->Uniform(-0.15, +0.15);
                //double vertex_y_smeared = gRandom->Uniform(-0.15, +0.15);
                //double vertex_z_smeared = z_target + gRandom->Uniform(-0.15, +0.15);
                //TVector3 vertex(vertex_x_smeared, vertex_y_smeared, vertex_z_smeared);
                TVector3 vertex(0, 0, z_target);
                TVector3 dir1 = (pos1 - vertex).Unit();
                TVector3 dir2 = (pos2 - vertex).Unit();

                TLorentzVector ph1(dir1.X() * ecal_e[icl], dir1.Y() * ecal_e[icl], dir1.Z() * ecal_e[icl], ecal_e[icl]);
                TLorentzVector ph2(dir2.X() * ecal_e[jcl], dir2.Y() * ecal_e[jcl], dir2.Z() * ecal_e[jcl], ecal_e[jcl]);


                Double_t pi0_mass = (ph1 + ph2).M();

                Double_t opening_angle = dir1.Angle(dir2) * (180.0 / TMath::Pi());
                if (opening_angle < 3.5 || opening_angle > 10) continue;  // The lower cut (e.g., < 6°) removes nearly collinear photon pairs → likely merged.
                                                                        // The upper cut (e.g., > 80°) removes highly unphysical, possibly misreconstructed pairs.


                if ((ph1 + ph2).M() <= 0.4 && (ph1 + ph2).M() >= 0){ 
                    // cuts to ensure some quality and avoid bugs (inf, -nan)
                    //cuts in time
                    if (clus_a_time[icl] > 100 && clus_a_time[icl] < 300 && 
                        clus_a_time[jcl] > 100 && clus_a_time[jcl] < 300 && 
                        fabs(clus_a_time[icl] - clus_a_time[jcl]) < 10) { // difference in time < 4 ns
                        h_pi0_mass->Fill(pi0_mass);
                        h_photons_time->Fill(clus_a_time[icl], clus_a_time[jcl]);
                        h_opening_angle->Fill(opening_angle);
                        h_xy->Fill(ecal_x[icl], ecal_y[icl]);
                        h_xy->Fill(ecal_x[jcl], ecal_y[jcl]);       
                    }
                        
                }
            }*/
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