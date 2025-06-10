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

void pi0_to_ntuple() {
    const Double_t z_calo = 6; // position of calorimeter from the target in m
    const Double_t z_target = 0.09;    // position of target
    //const Double_t z_origin = 0.0;
    //const Double_t vertex_z = z_target - z_origin;  // position of vertex, where the pi0 is created, right in the middle of the target

    TChain *ch = new TChain("T");

    //string path = "/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign_lowcur_April28/rootfiles/";
    //ifstream infile("lists_runs/runfiles_3173_to_3192.txt");

    //string path = "/adaqfs/home/a-onl/sbs/Rootfiles/";
    //ifstream infile("lists_runs/runfiles_3637.txt");

    // For the Monte Carlo simulation data
    string path = "/volatile/halla/sbs/efuchey/GEP3mod_pi0ECal_20250605/";
    ifstream infile("lists_runs/replays_mc.txt");

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

    Double_t goodblock_e[1000];
    Double_t goodblock_id[10000];
    Double_t goodblock_col[10000];
    Double_t goodblock_row[10000];
    Double_t goodblock_cid[10000];
    ch->SetBranchAddress("earm.ecal.goodblock.e", &goodblock_e);
    ch->SetBranchAddress("earm.ecal.goodblock.id", &goodblock_id);
    ch->SetBranchAddress("earm.ecal.goodblock.col", &goodblock_col);
    ch->SetBranchAddress("earm.ecal.goodblock.row", &goodblock_row);
    ch->SetBranchAddress("earm.ecal.goodblock.cid", &goodblock_cid);
    Double_t ngoodADChits = 0;
    ch->SetBranchAddress("earm.ecal.ngoodADChits", &ngoodADChits);

    TH1F *h_pi0_mass = new TH1F("h_pi0_mass", "Uncorrected #pi^{0} Invariant Mass;M_{#pi^{0}} [GeV];Events", 100, 0, 0.2);
    TH2F *h_photons_time = new TH2F("h_photons_time", "Time of photons;earm.ecal.clus.atimeblk[icl];earm.ecal.clus.atimeblk[jcl]", 100, 0, 540, 100, 0, 540);
    TH2F *h_xy = new TH2F("h_xy", "xy;earm.ecal.clus.x;earm.ecal.clus.y", 25, -1.5, 1.5, 25, -0.65, 0.6);
    TH1F *h_opening_angle = new TH1F("h_opening_angle", "Opening angle between #gamma#gamma;Angle [deg];Events", 100, 0, 30);

    // Output tree
    Double_t pi0_mass, opening_angle;
    Double_t ph1_t, ph2_t, ph1_x, ph1_y, ph2_x, ph2_y, ph1_E, ph2_E, ph1_px, ph1_py, ph1_pz, ph1_p, ph2_px, ph2_py, ph2_pz, ph2_p;
    Int_t n_blks = 0;
    Double_t nclus_out; // for output (avoid confusion with input nclus)
    Double_t blk_e[1000], blk_id[10000]; 
   
    TTree *output_tree = new TTree("pi0_tree", "Filtered #pi^{0} events");
    output_tree->Branch("invM", &pi0_mass, "invM/D");
    output_tree->Branch("ph1_t", &ph1_t, "ph1_t/D");
    output_tree->Branch("ph2_t", &ph2_t, "ph2_t/D");
    output_tree->Branch("ph1_x", &ph1_x, "ph1_x/D");
    output_tree->Branch("ph1_y", &ph1_y, "ph1_y/D");  
    output_tree->Branch("ph2_x", &ph2_x, "ph2_x/D");
    output_tree->Branch("ph2_y", &ph2_y, "ph2_y/D");
    output_tree->Branch("ph1_E", &ph1_E, "ph1_E/D");
    output_tree->Branch("ph2_E", &ph2_E, "ph2_E/D");
    output_tree->Branch("ph1_px", &ph1_px, "ph1_px/D");
    output_tree->Branch("ph1_py", &ph1_py, "ph1_py/D");
    output_tree->Branch("ph1_pz", &ph1_pz, "ph1_pz/D");
    output_tree->Branch("ph1_p", &ph1_p, "ph1_p/D");
    output_tree->Branch("ph2_px", &ph2_px, "ph2_px/D");
    output_tree->Branch("ph2_py", &ph2_py, "ph2_py/D");
    output_tree->Branch("ph2_pz", &ph2_pz, "ph2_pz/D");
    output_tree->Branch("ph2_p", &ph2_p, "ph2_p/D");
    output_tree->Branch("op_ang", &opening_angle, "op_ang/D");
    output_tree->Branch("nclus", &nclus_out, "nclus/D");
    output_tree->Branch("n_blks", &n_blks, "n_blks/I");
    output_tree->Branch("blk_e", &blk_e, "blk_e[1000]/D");
    output_tree->Branch("blk_id", &blk_id, "blk_id[10000]/D");


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
                //if (clus_a_time[icl] < 100 || clus_a_time[icl] > 140) continue;
                //if (clus_a_time[jcl] < 100 || clus_a_time[jcl] > 140) continue;
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
            ph1_t = clus_a_time[best_icl];
            ph2_t = clus_a_time[best_jcl];
            ph1_x = ecal_x[best_icl];
            ph1_y = ecal_y[best_icl];
            ph2_x = ecal_x[best_jcl];
            ph2_y = ecal_y[best_jcl];
            nclus_out = nclus;
            n_blks = clus_nblk[best_icl] + clus_nblk[best_jcl]; // total number of blocks in the two clusters
           
            for (int iblk = 0; iblk < ngoodADChits; ++iblk) {
                int cid = static_cast<int>(goodblock_cid[iblk]);
                // Check if the block belongs to one of the two clusters
                if (cid != best_icl && cid != best_jcl) continue; // skip blocks not in the selected clusters
                blk_e[iblk] = goodblock_e[iblk]; // copy the block energies to the output tree
                blk_id[iblk] = goodblock_id[iblk]; // copy the block IDs to the output tree

            }

            TVector3 pos1(ecal_x[best_icl], ecal_y[best_icl], z_calo);  // in m
            TVector3 pos2(ecal_x[best_jcl], ecal_y[best_jcl], z_calo);  // in m

            TVector3 vertex(0, 0, z_target);
            TVector3 dir1 = (pos1 - vertex).Unit();
            TVector3 dir2 = (pos2 - vertex).Unit();

            TLorentzVector ph1(dir1.X() * ecal_e[best_icl], dir1.Y() * ecal_e[best_icl], dir1.Z() * ecal_e[best_icl], ecal_e[best_icl]);
            TLorentzVector ph2(dir2.X() * ecal_e[best_jcl], dir2.Y() * ecal_e[best_jcl], dir2.Z() * ecal_e[best_jcl], ecal_e[best_jcl]);

            ph1_E = ecal_e[best_icl];
            ph2_E = ecal_e[best_jcl];
            ph1_px = dir1.X() * ecal_e[best_icl];
            ph1_py = dir1.Y() * ecal_e[best_icl];
            ph1_pz = dir1.Z() * ecal_e[best_icl];
            ph2_px = dir2.X() * ecal_e[best_jcl];
            ph2_py = dir2.Y() * ecal_e[best_jcl];
            ph2_pz = dir2.Z() * ecal_e[best_jcl];
            ph1_p = ph1.P(); // magnitude of the momentum vector
            ph2_p = ph2.P(); // magnitude of the momentum vector

            pi0_mass = (ph1 + ph2).M();

            opening_angle = dir1.Angle(dir2) * (180.0 / TMath::Pi());
            //if (opening_angle < 3) continue;
            if (opening_angle < 3.5 || opening_angle > 8) continue;  // The lower cut (e.g., < 6°) removes nearly collinear photon pairs → likely merged.
                                                                    // The upper cut (e.g., > 80°) removes highly unphysical, possibly misreconstructed pairs.
            //if (pi0_mass < 0 || pi0_mass > 0.4) continue;
            h_pi0_mass->Fill(pi0_mass);
            h_photons_time->Fill(clus_a_time[best_icl], clus_a_time[best_jcl]);
            h_opening_angle->Fill(opening_angle);
            h_xy->Fill(ecal_x[best_icl], ecal_y[best_icl]);
            h_xy->Fill(ecal_x[best_jcl], ecal_y[best_jcl]);
            // Fill the output tree
            output_tree->Fill();
        
        }    
        
    }

    TCanvas *c_pi0InvM = new TCanvas("c_pi0InvM", "Uncorrected #pi^{0} Invariant Mass", 800, 600);
    //h_pi0_mass->SetStats(0);
    h_pi0_mass->GetXaxis()->SetTitle("M_{#pi^{0}} [GeV]");
    h_pi0_mass->GetYaxis()->SetTitle("Events");
    h_pi0_mass->SetLineColor(kBlue);
    h_pi0_mass->SetMarkerStyle(20);
    h_pi0_mass->Draw();
    //c_pi0InvM->SaveAs("plots/pi0_invariant_mass_before.png");

    TCanvas *c_photons_time = new TCanvas("c_photons_time", "Time window", 800, 600);
    h_photons_time->SetStats(0);
    h_photons_time->GetXaxis()->SetTitle("t [ns]");
    h_photons_time->GetYaxis()->SetTitle("t [ns]");
    h_photons_time->SetTitle("Time of Photons");
    h_photons_time->Draw("colz");
    //c_photons_time->SaveAs("plots/photons_time.png");

    TCanvas *c_xy = new TCanvas("c_xy", "c_xy", 800, 600);
    h_xy->SetStats(0);
    h_xy->Draw("colz");
    h_xy->GetXaxis()->SetTitle("x [m]");
    h_xy->GetYaxis()->SetTitle("y [m]");
    //c_xy->SaveAs("plots/xy_clusters.png");

    TCanvas *c_opening = new TCanvas("c_opening", "Opening Angle", 800, 600);
    h_opening_angle->SetStats(0);
    h_opening_angle->GetXaxis()->SetTitle("Angle [deg]");
    h_opening_angle->GetYaxis()->SetTitle("Events");
    h_opening_angle->Draw();
    c_opening->SaveAs("plots/opening_angle.png");

    TFile *f = new TFile("histos/pi0_mass.root", "RECREATE");
    h_pi0_mass->Write();

    // Write the output tree to a file
    TFile *f_out = new TFile("histos/pi0_ntuple.root", "RECREATE");
    output_tree->Write();
    f_out->Close();
    cout << "Output tree with filtered #pi^{0} events written to histos/pi0_ntuple.root" << endl;
    cout << "Uncorrected #pi^{0} invariant mass histogram written to histos/pi0_mass.root" << endl;
    f->Close();
}