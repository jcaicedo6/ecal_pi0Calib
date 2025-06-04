#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TF1.h"
#include "TLegend.h"
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
#include "TDecompSVD.h"
//include random number generator
#include <TRandom3.h>



using namespace std;

// 
// run the script with "root -l ecal_pi0calib.C"

// Find the highest and second-highest energy clusters
void find_highest_energy(const Double_t ecal_e[], int nclus, int &imax1, int &imax2, double &e1, double &e2) {
    for (int icl = 0; icl < nclus; ++icl) {
        if (ecal_e[icl] > e1) {
            // If this cluster has the highest energy so far:
            e2 = e1; imax2 = imax1;               // The old highest becomes the second-highest, The old highest index becomes the second-highest index
            e1 = ecal_e[icl]; imax1 = icl;        // Update highest energy and highest index
        } else if (ecal_e[icl] > e2) {
            // If this cluster is not the highest, but is higher than the second-highest:
            e2 = ecal_e[icl]; imax2 = icl;        // Update second-highest energy and the second-highest index
        }
    }
}
// Best pair of clusters within the time window
void best_pair_cuts(
    const Double_t ecal_e[], const Double_t clus_nblk[],
    const Double_t ecal_x[], const Double_t ecal_y[], 
    const Double_t clus_a_time[], 
    int nclus, int &best_icl, int &best_jcl, double &best_dt) {
    for (int icl = 0; icl < nclus; ++icl) {
        for (int jcl = icl + 1; jcl < nclus; ++jcl) {
            if ((ecal_e[icl] + ecal_e[jcl]) < 0.5) continue;
            if (ecal_e[icl] < 0.2 || ecal_e[jcl] < 0.2) continue;
            if (clus_nblk[icl] < 2 || clus_nblk[jcl] < 2) continue;
            Double_t deltaR = sqrt(pow(ecal_x[icl] - ecal_x[jcl], 2) + pow(ecal_y[icl] - ecal_y[jcl], 2));
            if (deltaR < 0.07) continue;
            //pass_deltar++;
            if (clus_a_time[icl] < 100 || clus_a_time[icl] > 200) continue;
            if (clus_a_time[jcl] < 100 || clus_a_time[jcl] > 200) continue;
            double dt = fabs(clus_a_time[icl] - clus_a_time[jcl]);
            if (dt < best_dt && dt < 4) { // time window 
                best_dt = dt;
                best_icl = icl;
                best_jcl = jcl;
            }
        }
    }
}

void ecal_pi0calib() {
    const Double_t z_calo = 6; // position of calorimeter from the target in m
    const Double_t z_target = 0.09;    // position of target
    const Double_t z_origin = 0.0;
    const Double_t vertex_z = z_target - z_origin;  // position of vertex, where the pi0 is created, right in the middle of the target

    const int nbclusmax=100;// maximum number of clusters
    const int sizemax=1000;// maximum size of the cluster (how many blocks in each cluster)

    const double pi0_mass_pdg = 0.1349766;  // PDG pi0 mass in GeV
    TChain *ch = new TChain("T");

    //string path = "/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign_lowcur_April28/rootfiles/";
    string path = "/adaqfs/home/a-onl/sbs/Rootfiles/";
    //ifstream infile("lists_runs/runfiles_3173_to_3192.txt");
    ifstream infile("lists_runs/runfiles_3637.txt");
    //ifstream infile("lists_runs/runfiles_ecal.txt");
    string filename;

    while (getline(infile, filename)) {
        ch->Add((path + filename).c_str());
    }

    // Initialize the pointers before setting the branch address
    // Define the  energy, x, y of each cluster
    Double_t ecal_e[nbclusmax]; 
    Double_t ecal_x[nbclusmax];
    Double_t ecal_y[nbclusmax];

    ch->SetBranchStatus("*", 0); // disable all Branches
    ch->SetBranchAddress("earm.ecal.clus.e", &ecal_e);
    ch->SetBranchAddress("earm.ecal.clus.x", &ecal_x);
    ch->SetBranchAddress("earm.ecal.clus.y", &ecal_y);

    // Define the number of clusters
    Double_t nclus = 0;
    ch->SetBranchAddress("earm.ecal.nclus", &nclus);

    // Define the time of each cluster
    Double_t clus_a_time[1000];
    ch->SetBranchAddress("earm.ecal.clus.atimeblk", &clus_a_time);

    // Define the number of blocks in each cluster and the id of each block per cluster (max 100 blocks per cluster)
    Double_t clus_nblk[nbclusmax];       // number of blocks in the cluster
    Double_t clus_id[nbclusmax];         // block id
    Double_t clus_row[nbclusmax];        // block row in the calorimeter
    Double_t clus_col[nbclusmax];        // block column in the calorimeter
    //Double_t clus_eblk[nbclusmax];       // block energy with highest energy in the cluster
    //Double_t clus_idblk[nbclusmax];      // block id with highest energy in the cluster
    Double_t goodblock_e[1000];
    Double_t goodblock_id[10000];
    Double_t goodblock_col[10000];
    Double_t goodblock_row[10000];
    Double_t goodblock_cid[10000];
    //Double_t energy_blk[1000];

    ch->SetBranchAddress("earm.ecal.clus.nblk", &clus_nblk);
    ch->SetBranchAddress("earm.ecal.clus.id", &clus_id);
    ch->SetBranchAddress("earm.ecal.clus.row", &clus_row);
    ch->SetBranchAddress("earm.ecal.clus.col", &clus_col);
    //ch->SetBranchAddress("earm.ecal.clus.eblk", &clus_eblk);
    //ch->SetBranchAddress("earm.ecal.idblk", &clus_idblk);
    ch->SetBranchAddress("earm.ecal.goodblock.e", &goodblock_e);
    ch->SetBranchAddress("earm.ecal.goodblock.id", &goodblock_id);
    ch->SetBranchAddress("earm.ecal.goodblock.col", &goodblock_col);
    ch->SetBranchAddress("earm.ecal.goodblock.row", &goodblock_row);
    ch->SetBranchAddress("earm.ecal.goodblock.cid", &goodblock_cid);


    Double_t ngoodADChits = 0;
    ch->SetBranchAddress("earm.ecal.ngoodADChits", &ngoodADChits);

    TH1F *h_pi0_mass = new TH1F("h_pi0_mass", "Pi0 Invariant Mass;M_{#pi^{0}} [GeV];Events", 80, 0, 0.6);
    TH1F *h_pi0_mass_corr = new TH1F("h_pi0_mass_reco", "Reconstructed Pi0 Invariant Mass;M_{#pi^{0}} [GeV];Events", 80, 0, 0.6);

    Long64_t nEvents = ch->GetEntries();
    cout << "Number of events: " << nEvents << endl;


    // --- Discover geometry from data: collect all unique block IDs ---
    int minRow = INT_MAX, maxRow = INT_MIN;
    int minCol = INT_MAX, maxCol = INT_MIN;
    std::set<int> blockIDs;
    std::map<int, int> blockID_to_row, blockID_to_col;

    for (Long64_t i = 0; i < nEvents; ++i) {
        ch->GetEntry(i);
        if (nclus < 1) continue;
        for (unsigned int b = 0; b < (unsigned int)ngoodADChits; ++b) {
            int bid = static_cast<int>(goodblock_id[b]);
            int row = static_cast<int>(goodblock_row[b]);
            int col = static_cast<int>(goodblock_col[b]);
            if (bid >= 0) {
                blockIDs.insert(bid);
                blockID_to_row[bid] = row;
                blockID_to_col[bid] = col;
                if (row < minRow) minRow = row;
                if (row > maxRow) maxRow = row;
                if (col < minCol) minCol = col;
                if (col > maxCol) maxCol = col;
            }
        }
    }
    int nlin = maxRow - minRow + 1;
    int ncol = maxCol - minCol + 1;

    for (Long64_t i = 0; i < nEvents/100; ++i) {
        ch->GetEntry(i);
        if (nclus < 1) continue;
        // All good blocks (only active blocks)
        for (unsigned int b = 0; b < (unsigned int)ngoodADChits; ++b) {
            int bid = static_cast<int>(goodblock_id[b]);
            if (bid >= 0) blockIDs.insert(bid);
        }
    }

    // Build mapping: blockID <-> matrix index
    std::map<int, int> blockID_to_idx;
    std::vector<int> idx_to_blockID;
    int idx = 0;
    for (int bid : blockIDs) {
        blockID_to_idx[bid] = idx++;
        idx_to_blockID.push_back(bid);
    }
    int nblocks = idx_to_blockID.size();

    cout << "Number of unique blocks: " << nblocks << endl;

    // Matrix used for calibration
    TMatrixD A(nblocks, nblocks);
    TVectorD B(nblocks);
    A.Zero();
    B.Zero();
  
    // Block occupancy for later use
    std::vector<int> occupancy(nblocks, 0);

    // Define the counters for the cuts
    Double_t pass_clus=0, pass_deltar=0, pass_time=0, pass_mass=0;
    // initialize the random number generator
    gRandom->SetSeed(0);
    //Double_t pi0_mass_pdg_random = pi0_mass_pdg;    
    std::vector<std::vector<std::vector<int>>> clusterBlocks(nEvents);

    for (Long64_t i = 0; i < nEvents; i++) {
        ch->GetEntry(i);
        if (i % 10000 == 0) {
            cout << "Processing event (correction) " << i << " / " << nEvents << "\r";
            cout.flush();
        }
        // Check that there are at least two clusters
        if (nclus < 2) continue;
        pass_clus++;
        // Find the two clusters with the best time difference and apply cuts
        double best_dt = 1e9;
        int best_icl = -1, best_jcl = -1;
        best_pair_cuts(ecal_e, clus_nblk, ecal_x, ecal_y, clus_a_time, nclus, best_icl, best_jcl, best_dt);
        pass_time++;
        // Find indices of the two clusters with the highest energies
        int imax1 = -1, imax2 = -1; // Index of highest and second-highest energy cluster
        double e1 = -1, e2 = -1;  // Energy of highest and second-highest energy cluster
        find_highest_energy(ecal_e, nclus, imax1, imax2, e1, e2);

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
            if (pi0_mass < 0.02 || pi0_mass > 0.4) continue;
            pass_mass++;
            h_pi0_mass->Fill(pi0_mass);

            // compute expected total π0 energy as using the scale factor
            //double pi0_mass_smeared = gRandom->Gaus(pi0_mass_pdg, 0.005); // 5 MeV width max
            double expected_E = (ecal_e[best_icl] + ecal_e[best_jcl]) * (pi0_mass_pdg / pi0_mass);
            //double expected_E = pi0_mass_pdg / pi0_mass;
            
            // **Loop over all blocks in cluster **:
            // zero the per-block energy accumulator
            std::vector<double> energy(nblocks, 0.0);
    
            // Loop over all clusters and blocks per event
            for (int cl = 0; cl < nclus; ++cl) {
                // Assign energy to each good block
                for (unsigned int b = 0; b < (unsigned int)ngoodADChits; ++b) {
                    int rawID = (int)goodblock_id[b];
                    int cid = (int)goodblock_cid[b];
                    //int iblock = num_btom[rawID];
                    int iblock = blockID_to_idx.count(rawID) ? blockID_to_idx[rawID] : -1;
                    if (iblock < 0 || iblock >= nblocks) continue;
                    if (cid != cl) continue; // Only assign energy to blocks in this cluster
                    energy[iblock] += goodblock_e[b];
                    if (energy[iblock] > 0.0) occupancy[iblock]++;
                }
            }
            // Now fill A and B:
            for (int j = 0; j < nblocks; ++j) {
                for (int k = 0; k < nblocks; ++k) {
                    A(j,k) += energy[j] * energy[k];
                }
                B(j) += energy[j] * expected_E;
            } 
        }
    }
    cout<<"Events passing number of clusters cut: "<<pass_clus<<endl;
    //cout<<"Events passing deltaR cut: "<<pass_deltar<<endl;
    cout<<"Events passing time cut: "<<pass_time<<endl;
    cout<<"Events passing all cuts: "<<pass_mass<<endl;
  
    // ===================== OCCUPANCY PRUNING & REGULARIZATION =====================
    const int N_min = 30;
    for (int i = 0; i < nblocks; ++i) {
        if (occupancy[i] < N_min) {
            // zero out row and column i
            for (int j = 0; j < nblocks; ++j) {
                A(i,j) = A(j,i) = 0.0;
            }
            //B(0,i) = 0.0;
            B(i) = 0.0;
        }
    }
    // add small diagonal term λ·A_ii for singularity
    const double lambda = 1e-4;
    for (int i = 0; i < nblocks; ++i) {
        A(i,i) += lambda * A(i,i);
    }
        
    // -------------------- Invert the Matrix --------------------
    //A.Print();
    std::cout
        <<"[DEBUG] A(0,0)="<<A(0,0)
        <<", B(0,0)="<<B(0)
        <<"\n";
        for (int i = 0; i < std::min(5,nblocks); ++i) {
        std::cout<<"  B(0,"<<i<<")="<<B(i)<<"\n";
        }

    //Double_t coeff[nblocks];
    std::vector<double> coeff(nblocks, 1.0);

    // SVD decomposition of A and solve C = A⁺·B (see https://root.cern.ch/doc/master/classTDecompSVD.html)
    TDecompSVD svd(A);
    if (!svd.Decompose()) {
      cerr<<"ERROR: SVD decomposition failed\n";
      return;
    }
    // Prepare solution vector
    // Copy B into C, then solve C = A⁺·B
    TVectorD C = B;           // initialize C with RHS
    Bool_t ok = svd.Solve(C); // Solve in place: C ← pseudo-inverse(A)·C
    if (!ok) {
      cerr<<"ERROR: SVD Solve failed\n";
      return;
    }

    for (int i = 0; i < min(nblocks,5); ++i) {
      cout << "C["<<i<<"] = " << C(i) << "\n";
    }
        // the coefficients of the linear system 
        // C = A^-1 * B where A^-1 is the inverse of A, B is the vector of unknowns, and C is the vector of calibration coefficients. originally A * C = B
        //TMatrixD C = B * A;  // 1×nblocksm * nblocksm×nblocksm

    // now overwrite just the “good” ones using occupancy
    for (int compID = 0; compID < nblocks; ++compID) {
        //int rawID = num_mtob[compID];     // reverse map
        //int blockID = idx_to_blockID[compID];
        //double c = C(0,compID);
        double c = C(compID);

        if (occupancy[compID] < N_min || c < 1e-2 || c > 10) {
            coeff[compID] = 1.0;
        } else {
            coeff[compID] = c;
        //coeff[rawID] = C(0, compID);      // use same compID as above
        }
    }
    
    // -------------------- Saving the coefficients in a text file --------------------
    cout<<"Writing coefficients to file..."<<endl;

    ofstream coeff_file("ecal_block_calibration_factors.txt");
    coeff_file << "# BlockID\tCalibrationCoeff\n";

    for (int i = 0; i < nblocks; i++) {
        int blockID = idx_to_blockID[i];
        coeff_file << blockID << "\t" << coeff[i] << endl;
        cout << "BlockID: " << i << "\tCalibrationCoeff: " << coeff[i] << endl;
    }

    coeff_file.close();

    // Now, with the coefficients in hand, I can use them to calibrate the Ecal energy in the code.
    // Use the coefficients to calibrate the Ecal energy and plot the pi0 invariant mass
    // with and without the correction

    // counters for passing events
    Double_t pass_clus_corr=0, pass_deltar_corr=0, pass_time_corr=0, pass_mass_corr=0;
    double test_mass_sum = 0;
    int test_mass_count = 0;
 
    for (Long64_t i = 0; i < nEvents; i++) {
        ch->GetEntry(i);
        if (i % 10000 == 0) {
            cout << "Processing event (correction) " << i << " / " << nEvents << "\r";
            cout.flush();
        }
        if (nclus < 2) continue;
        //pass_clus_corr++;
         // Find the two clusters with the best time difference and apply cuts
        double best_dt = 1e9;
        int best_icl = -1, best_jcl = -1;
        best_pair_cuts(ecal_e, clus_nblk, ecal_x, ecal_y, clus_a_time, nclus, best_icl, best_jcl, best_dt);
        //pass_time_corr++;
        // Find indices of the two clusters with the highest energies
        int imax1 = -1, imax2 = -1; // Index of highest and second-highest energy cluster
        double e1 = -1, e2 = -1;  // Energy of highest and second-highest energy cluster
        find_highest_energy(ecal_e, nclus, imax1, imax2, e1, e2);
 
        // **After** correction: rebuild each photon’s energy by summing
        //     per‑block energies * coefficients
        if (best_icl >= 0 && best_jcl >= 0) {
            if (!((best_icl == imax1 && best_jcl == imax2) || (best_icl == imax2 && best_jcl == imax1))) {
                continue; // skip if not the two highest-energy clusters
            }
            double e_corr[2] = {0., 0.};

            int idx[2] = {best_icl, best_jcl};
            for (int ic = 0; ic < 2; ++ic) {
                int cl = idx[ic];

                // Assign energy to each good block
                for (unsigned int b = 0; b < (unsigned int)ngoodADChits; ++b) {
                    int rawID = (int)goodblock_id[b];
                    int cid = (int)goodblock_cid[b];
                    // Use blockID_to_idx to get the matrix index
                    auto it = blockID_to_idx.find(rawID);
                    if (it == blockID_to_idx.end()) continue; // skip if not a known block
                    if (cid != cl) continue; // Only assign energy to blocks in this cluster
                    e_corr[ic] += coeff[it->second] * goodblock_e[b];
                }
            }
            TVector3 pos1(ecal_x[best_icl], ecal_y[best_icl], z_calo);  // in m
            TVector3 pos2(ecal_x[best_jcl], ecal_y[best_jcl], z_calo);  // in m

            TVector3 vertex(0, 0, z_target);
            TVector3 dir1 = (pos1 - vertex).Unit();
            TVector3 dir2 = (pos2 - vertex).Unit();

            // rebuild corrected TLorentzVectors
            TLorentzVector ph1_corr(dir1.X() * e_corr[0], dir1.Y() * e_corr[0], dir1.Z() * e_corr[0], e_corr[0]);
            TLorentzVector ph2_corr(dir2.X() * e_corr[1], dir2.Y() * e_corr[1], dir2.Z() * e_corr[1], e_corr[1]);
            
            double pi0_mass_corr = (ph1_corr + ph2_corr).M();

            Double_t opening_angle = dir1.Angle(dir2) * (180.0 / TMath::Pi());
            if (opening_angle < 3.5 || opening_angle > 8) continue;  // The lower cut (e.g., < 6°) removes nearly collinear photon pairs → likely merged.
                                                                    // The upper cut (e.g., > 80°) removes highly unphysical, possibly misreconstructed pairs.

            if (pi0_mass_corr <= 0.02 || pi0_mass_corr >= 0.4) continue;    // Making a cut on the pi0 mass, e.g., between 0.06 and 0.6 GeV

            // fill histograms
            test_mass_sum += pi0_mass_corr;
            test_mass_count++;
            //pass_mass_corr++;
            //h_pi0_mass_corr->Fill(pi0_mass_corr);
        }
    }
    cout << "Scaling coefficients to match PDG π⁰ mass..." << endl;
    if (test_mass_sum > 0) {
        double mean_mass = test_mass_sum / test_mass_count;
        double scale_factor = pi0_mass_pdg / mean_mass;

        // Apply the scaling factor to the coefficients
        for (int i = 0; i < nblocks; ++i) {
            coeff[i] *= scale_factor;
        }
    }
    //h_pi0_mass_corr->Reset();

    for (Long64_t i = 0; i < nEvents; i++) {
        ch->GetEntry(i);
        cout << "Processing event " << i << "\r";
        if (nclus < 2) continue;
        pass_clus_corr++;
         // Find the two clusters with the best time difference and apply cuts
        double best_dt = 1e9;
        int best_icl = -1, best_jcl = -1;
        best_pair_cuts(ecal_e, clus_nblk, ecal_x, ecal_y, clus_a_time, nclus, best_icl, best_jcl, best_dt);
        pass_time_corr++;
        // Find indices of the two clusters with the highest energies
        int imax1 = -1, imax2 = -1; // Index of highest and second-highest energy cluster
        double e1 = -1, e2 = -1;  // Energy of highest and second-highest energy cluster
        find_highest_energy(ecal_e, nclus, imax1, imax2, e1, e2);
 
        // **After** correction: rebuild each photon’s energy by summing
        //     per‑block energies * coefficients
        if (best_icl >= 0 && best_jcl >= 0) {
            if (!((best_icl == imax1 && best_jcl == imax2) || (best_icl == imax2 && best_jcl == imax1))) {
                continue; // skip if not the two highest-energy clusters
            }
            double e_corr[2] = {0., 0.};

            int idx[2] = {best_icl, best_jcl};
            for (int ic = 0; ic < 2; ++ic) {
                int cl = idx[ic];

                // Assign energy to each good block
                for (unsigned int b = 0; b < (unsigned int)ngoodADChits; ++b) {
                    int rawID = (int)goodblock_id[b];
                    int cid = (int)goodblock_cid[b];
                    // Use blockID_to_idx to get the matrix index
                    auto it = blockID_to_idx.find(rawID);
                    if (it == blockID_to_idx.end()) continue; // skip if not a known block
                    if (cid != cl) continue; // Only assign energy to blocks in this cluster
                    e_corr[ic] += coeff[it->second] * goodblock_e[b];
                }
            }
            TVector3 pos1(ecal_x[best_icl], ecal_y[best_icl], z_calo);  // in m
            TVector3 pos2(ecal_x[best_jcl], ecal_y[best_jcl], z_calo);  // in m

            TVector3 vertex(0, 0, z_target);
            TVector3 dir1 = (pos1 - vertex).Unit();
            TVector3 dir2 = (pos2 - vertex).Unit();

            // rebuild corrected TLorentzVectors
            TLorentzVector ph1_corr(dir1.X() * e_corr[0], dir1.Y() * e_corr[0], dir1.Z() * e_corr[0], e_corr[0]);
            TLorentzVector ph2_corr(dir2.X() * e_corr[1], dir2.Y() * e_corr[1], dir2.Z() * e_corr[1], e_corr[1]);
            
            double pi0_mass_corr = (ph1_corr + ph2_corr).M();

            Double_t opening_angle = dir1.Angle(dir2) * (180.0 / TMath::Pi());
            if (opening_angle < 3.5 || opening_angle > 8) continue;  // The lower cut (e.g., < 6°) removes nearly collinear photon pairs → likely merged.
                                                                    // The upper cut (e.g., > 80°) removes highly unphysical, possibly misreconstructed pairs.

            if (pi0_mass_corr <= 0.02 || pi0_mass_corr >= 0.4) continue;    // Making a cut on the pi0 mass, e.g., between 0.06 and 0.6 GeV

            // fill histograms
            test_mass_sum += pi0_mass_corr;
            test_mass_count++;
            pass_mass_corr++;
            h_pi0_mass_corr->Fill(pi0_mass_corr);
        }
    }

    // -------------------- Plotting --------------------
    TCanvas *c = new TCanvas("c","before/after", 800,600);
    h_pi0_mass->SetLineColor(kRed);
    h_pi0_mass_corr->SetLineColor(kBlue);
    h_pi0_mass_corr->Draw();
    h_pi0_mass->Draw("SAME");
    auto leg = new TLegend(0.6,0.7,0.9,0.9);
    leg->AddEntry(h_pi0_mass,      "Before calib", "l");
    leg->AddEntry(h_pi0_mass_corr, "After calib",  "l");
    leg->Draw();


    // Create a 2D histogram to visualize calibration coefficients in the detector geometry
    TH2F *h_coeff_map = new TH2F("h_coeff_map", "ECAL Block Calibration Coefficients;Column;Row",
        ncol, minCol, maxCol + 1l,
        nlin, minRow, maxRow + 1);

    // Fill the 2D histogram using the coeff[] array 
    for (int i = 0; i < nblocks; ++i) {
        int blockID = idx_to_blockID[i];
        int row = blockID_to_row[blockID];
        int col = blockID_to_col[blockID];
        double val = coeff[i];
        h_coeff_map->SetBinContent(col - minCol + 1, row - minRow + 1, val);  // ROOT bins start from 1
    }
    cout<<"Events passing number of clusters cut after calib: "<<pass_clus_corr<<endl;
    //cout<<"Events passing deltaR cut after calib: "<<pass_deltar_corr<<endl;
    cout<<"Events passing time cut after calib: "<<pass_time_corr<<endl;
    cout<<"Events passing all cuts after calib: "<<pass_mass_corr<<endl;

    cout<<"-------------------- Difference between before and after calib -------------------"<<endl;
    cout << "Before calib: mean = " << h_pi0_mass->GetMean() 
        << ", sigma = " << h_pi0_mass->GetRMS() << endl;
    cout << "After calib:  mean = " << h_pi0_mass_corr->GetMean() 
        << ", sigma = " << h_pi0_mass_corr->GetRMS() << endl;

    // Draw the map
    TCanvas *c_map = new TCanvas("c_map", "ECAL Coefficients Heatmap", 900, 800);
    h_coeff_map->SetStats(0);
    h_coeff_map->Draw("COLZ");
}