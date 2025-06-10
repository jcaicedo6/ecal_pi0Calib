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
#include <vector>
#include <string>
#include "TCanvas.h"
// for ifstream
#include <fstream>
#include "TDecompSVD.h"
//include random number generator
#include <TRandom3.h>
#include <set>
#include <map>
#include <algorithm>

using namespace std;

// Find the highest and second-highest energy clusters
void find_highest_energy(const Double_t ecal_e[], int nclus, int &imax1, int &imax2, double &e1, double &e2) {
    imax1 = -1; imax2 = -1;
    e1 = -1; e2 = -1;
    for (int icl = 0; icl < nclus; ++icl) {
        if (ecal_e[icl] > e1) {
            e2 = e1; imax2 = imax1;
            e1 = ecal_e[icl]; imax1 = icl;
        } else if (ecal_e[icl] > e2) {
            e2 = ecal_e[icl]; imax2 = icl;
        }
    }
}

// Best pair of clusters within the time window
void best_pair_cuts(
    const Double_t ecal_e[], const Double_t clus_nblk[],
    const Double_t ecal_x[], const Double_t ecal_y[],
    const Double_t clus_a_time[],
    int nclus, int &best_icl, int &best_jcl, double &best_dt) {
    best_icl = -1;
    best_jcl = -1;
    best_dt = 1e9;
    for (int icl = 0; icl < nclus; ++icl) {
        for (int jcl = icl + 1; jcl < nclus; ++jcl) {
            if ((ecal_e[icl] + ecal_e[jcl]) < 0.5) continue;
            if (ecal_e[icl] < 0.2 || ecal_e[jcl] < 0.2) continue;
            if (clus_nblk[icl] < 2 || clus_nblk[jcl] < 2) continue;
            Double_t deltaR = sqrt(pow(ecal_x[icl] - ecal_x[jcl], 2) + pow(ecal_y[icl] - ecal_y[jcl], 2));
            if (deltaR < 0.07) continue;
            if (clus_a_time[icl] < 100 || clus_a_time[icl] > 140) continue;
            if (clus_a_time[jcl] < 100 || clus_a_time[jcl] > 140) continue;
            double dt = fabs(clus_a_time[icl] - clus_a_time[jcl]);
            if (dt < best_dt && dt < 4) { // time window
                best_dt = dt;
                best_icl = icl;
                best_jcl = jcl;
            }
        }
    }
}

void ecal_pi0calib_iter(int n_iterations = 2) {
    const Double_t z_calo = 6;
    const Double_t z_target = 0.09;
    //const Double_t z_origin = 0.0;
    //const Double_t vertex_z = z_target - z_origin;
    TVector3 vertex(0, 0, z_target);

    const int nbclusmax = 100;
    const int sizemax = 1000;
    const double pi0_mass_pdg = 0.1349766;

    TChain *ch = new TChain("T");
    string path = "/adaqfs/home/a-onl/sbs/Rootfiles/";
    ifstream infile("lists_runs/runfiles_3637.txt");
    string filename;
    while (getline(infile, filename)) {
        ch->Add((path + filename).c_str());
    }

    Double_t ecal_e[nbclusmax];
    Double_t ecal_x[nbclusmax];
    Double_t ecal_y[nbclusmax];
    Double_t nclus = 0;
    Double_t clus_a_time[nbclusmax]; // Corrected size
    Double_t clus_nblk[nbclusmax];
    //Double_t clus_id[nbclusmax];
    //Double_t clus_row[nbclusmax];
    //Double_t clus_col[nbclusmax];
    Double_t goodblock_e[sizemax];
    Double_t goodblock_id[sizemax * 10];
    Double_t goodblock_col[sizemax * 10];
    Double_t goodblock_row[sizemax * 10];
    Double_t goodblock_cid[sizemax * 10];
    Double_t ngoodADChits = 0;

    ch->SetBranchStatus("*", 0);
    ch->SetBranchAddress("earm.ecal.clus.e", &ecal_e);
    ch->SetBranchAddress("earm.ecal.clus.x", &ecal_x);
    ch->SetBranchAddress("earm.ecal.clus.y", &ecal_y);
    ch->SetBranchAddress("earm.ecal.nclus", &nclus);
    ch->SetBranchAddress("earm.ecal.clus.atimeblk", &clus_a_time);
    ch->SetBranchAddress("earm.ecal.clus.nblk", &clus_nblk);
    //ch->SetBranchAddress("earm.ecal.clus.id", &clus_id);
    //ch->SetBranchAddress("earm.ecal.clus.row", &clus_row);
    //ch->SetBranchAddress("earm.ecal.clus.col", &clus_col);
    ch->SetBranchAddress("earm.ecal.goodblock.e", &goodblock_e);
    ch->SetBranchAddress("earm.ecal.goodblock.id", &goodblock_id);
    ch->SetBranchAddress("earm.ecal.goodblock.col", &goodblock_col);
    ch->SetBranchAddress("earm.ecal.goodblock.row", &goodblock_row);
    ch->SetBranchAddress("earm.ecal.goodblock.cid", &goodblock_cid);
    ch->SetBranchAddress("earm.ecal.ngoodADChits", &ngoodADChits);

    TH1F *h_pi0_mass = new TH1F("h_pi0_mass", "Reconstructed #pi^{0} Invariant Mass;M_{#pi^{0}} [GeV];Events", 80, 0, 0.6);
    TH1F *h_pi0_mass_corr = new TH1F("h_pi0_mass_reco", "Reconstructed #pi^{0} Invariant Mass;M_{#pi^{0}} [GeV];Events", 80, 0, 0.6);

    Long64_t nEvents = ch->GetEntries();
    cout << "Number of events: " << nEvents << endl;

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
            }
        }
    }
    int minRow = INT_MAX, maxRow = INT_MIN;
    int minCol = INT_MAX, maxCol = INT_MIN;
    for (const auto& pair : blockID_to_row) {
        minRow = min(minRow, pair.second);
        maxRow = max(maxRow, pair.second);
    }
    for (const auto& pair : blockID_to_col) {
        minCol = min(minCol, pair.second);
        maxCol = max(maxCol, pair.second);
    }
    int nlin = maxRow - minRow + 1;
    int ncol = maxCol - minCol + 1;

    std::map<int, int> blockID_to_idx;
    std::vector<int> idx_to_blockID;
    int idx = 0;
    for (int bid : blockIDs) {
        blockID_to_idx[bid] = idx++;
        idx_to_blockID.push_back(bid);
    }
    int nblocks = idx_to_blockID.size();
    cout << "Number of unique blocks: " << nblocks << endl;

    // This vector stores the calibration coefficients, updated iteratively.
    std::vector<double> current_coeff(nblocks, 1.0);

    // --- Iterative Calibration Loop ---
    for (int iter = 0; iter < n_iterations; ++iter) {
        cout << "\n--- Iteration " << iter + 1 << " ---" << endl;
        TMatrixD A(nblocks, nblocks);
        TVectorD B(nblocks);
        std::vector<int> occupancy(nblocks, 0);
        A.Zero();
        B.Zero();
        Double_t pass_clus = 0, pass_time = 0, pass_mass = 0;
        double sum_mass_uncorrected = 0;
        int count_mass_uncorrected = 0;

        for (Long64_t i = 0; i < nEvents; i++) {
            ch->GetEntry(i);
            if (i % 10000 == 0) {
                cout << "Processing event " << i << " / " << nEvents << "\r";
                cout.flush();
            }
            if (nclus < 2) continue;
            pass_clus++;

            double best_dt = 1e9;
            int best_icl = -1, best_jcl = -1;
            best_pair_cuts(ecal_e, clus_nblk, ecal_x, ecal_y, clus_a_time, nclus, best_icl, best_jcl, best_dt);

            int imax1 = -1, imax2 = -1;
            double e1 = -1, e2 = -1;
            find_highest_energy(ecal_e, nclus, imax1, imax2, e1, e2);

            if (best_icl >= 0 && best_jcl >= 0 && ((best_icl == imax1 && best_jcl == imax2) || (best_icl == imax2 && best_jcl == imax1))) {
                pass_time++;
                TVector3 pos1(ecal_x[best_icl], ecal_y[best_icl], z_calo);
                TVector3 pos2(ecal_x[best_jcl], ecal_y[best_jcl], z_calo);
                //TVector3 vertex(0, 0, z_target);
                TVector3 dir1 = (pos1 - vertex).Unit();
                TVector3 dir2 = (pos2 - vertex).Unit();
                TLorentzVector ph1(dir1.X() * ecal_e[best_icl], dir1.Y() * ecal_e[best_icl], dir1.Z() * ecal_e[best_icl], ecal_e[best_icl]);
                TLorentzVector ph2(dir2.X() * ecal_e[best_jcl], dir2.Y() * ecal_e[best_jcl], dir2.Z() * ecal_e[best_jcl], ecal_e[best_jcl]);
                Double_t pi0_mass = (ph1 + ph2).M();
                Double_t opening_angle = dir1.Angle(dir2) * (180.0 / TMath::Pi());

                if (opening_angle < 3.5 || opening_angle > 8) continue;
                if (pi0_mass > 0.02 && pi0_mass < 0.4) {
                    pass_mass++;
                    sum_mass_uncorrected += pi0_mass;
                    count_mass_uncorrected++;
                    if (iter == 0) {
                        // Fill the uncorrected mass histogram only in the first iteration
                        h_pi0_mass->Fill(pi0_mass);
                    } 

                    double expected_E = (ecal_e[best_icl] + ecal_e[best_jcl]) * (pi0_mass_pdg / pi0_mass);
                    static std::vector<double> energy(nblocks, 0.0); // Estimated energy per unique block for this event
                    energy.assign(nblocks, 0.0);    // Reset energy accumulator for this event

                    for (unsigned int b = 0; b < (unsigned int)ngoodADChits; ++b) {
                        int rawID = (int)goodblock_id[b];
                        int cid = (int)goodblock_cid[b];
                        // Map raw block ID to its index in the calibration matrix
                        int iblock = blockID_to_idx.count(rawID) ? blockID_to_idx[rawID] : -1;
                        if (iblock >= 0 && iblock < nblocks) {
                            bool in_best_clusters = false;
                            // Check if this goodblock belongs to one of the two best clusters
                            for (int cl_idx : {best_icl, best_jcl}) {
                                if (cid == cl_idx) {
                                    in_best_clusters = true;
                                    break;
                                }
                            }
                            if (in_best_clusters) {
                                // Apply the calibration coefficient from the previous iteration
                                energy[iblock] += goodblock_e[b];
                                if (energy[iblock] > 0.0) occupancy[iblock]++;
                            }
                        }
                    }

                    // Accumulate terms for the A matrix and B vector
                    for (int j = 0; j < nblocks; ++j) {
                        for (int k = 0; k < nblocks; ++k) {
                            A(j, k) += energy[j] * energy[k];
                        }
                        B(j) += energy[j] * expected_E;
                    }
                }
            }
        }
        cout << "Iteration " << iter + 1 << " - Events passing cuts: Clusters=" << pass_clus << ", Time=" << pass_time << ", Mass=" << pass_mass << endl;
        cout << "Iteration " << iter + 1 << " - Mean uncorrected mass: ";
        if (count_mass_uncorrected > 0) cout << sum_mass_uncorrected / count_mass_uncorrected << endl; else cout << "N/A" << endl;

        // ===================== OCCUPANCY PRUNING & REGULARIZATION =====================
        const int N_min = 30; // Minimum occupancy for a block to be calibrated
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
        const double lambda = 1e-4;  // Regularization parameter
        for (int i = 0; i < nblocks; ++i) {
            A(i,i) += lambda * A(i,i);
        }
        // --- Solve the Linear System (A * C = B) using SVD ---
        TDecompSVD svd(A);
        if (!svd.Decompose()) {
            cerr << "ERROR: SVD decomposition failed in iteration " << iter + 1 << endl;
            return;
        }
        TVectorD C = B; // Initialize solution vector with B
        Bool_t ok = svd.Solve(C); // Solve in place: C becomes pseudo-inverse(A) * B
        if (!ok) {
            cerr << "ERROR: SVD Solve failed in iteration " << iter + 1 << endl;
            // Continue, but coefficients might be problematic
        }

        // --- Update Coefficients for Next Iteration ---
        double sum_new_coeff = 0;
        int count_new_coeff = 0;
        double sum_delta = 0;
        double max_delta = 0;
        std::vector<double> new_coeff(nblocks, 1.0); // Temporary vector for new coefficients
        for (int i = 0; i < nblocks; ++i) {
            // Only update if block meets occupancy and coefficient is within reasonable bounds
            if (occupancy[i] >= N_min && C(i) > 1e-2 && C(i) < 10) {
                new_coeff[i] = C(i);
                sum_new_coeff += C(i);
                count_new_coeff++;
                double delta = fabs(C(i) - current_coeff[i]);
                sum_delta += delta;
                if (delta > max_delta) max_delta = delta;
            } else {
                new_coeff[i] = current_coeff[i]; // Keep previous coefficient if new one is bad
            }
        }
        cout << "Iteration " << iter + 1 << " - Average new coefficient: ";
        if (count_new_coeff > 0) cout << sum_new_coeff / count_new_coeff << endl; else cout << "N/A" << endl;

        cout << "Iteration " << iter + 1 << " - Average coefficient change: ";
        if (count_new_coeff > 0) cout << sum_delta / count_new_coeff << endl; else cout << "N/A" << endl;

        cout << "Iteration " << iter + 1 << " - Max coefficient change: " << max_delta << endl;

        current_coeff = new_coeff; // Update the main coefficients vector for the next iteration

        // Save coefficients for each iteration
        ofstream coeff_file("ecal_block_calibration_factors_iter_" + to_string(iter + 1) + ".txt");
        coeff_file << "# BlockID\tCalibrationCoeff\n";
        for (int i = 0; i < nblocks; i++) {
            int blockID = idx_to_blockID[i];
            coeff_file << blockID << "\t" << current_coeff[i] << endl;
        }
        coeff_file.close();
        cout << "Coefficients for iteration " << iter + 1 << " saved." << endl;
    }

    // --- Apply final calibration and plot ---
    // Reset the corrected histogram before filling it with the final calibrated data
    h_pi0_mass_corr->Reset();
    Double_t pass_clus_corr = 0, pass_time_corr = 0, pass_mass_corr = 0;
    double test_mass_sum = 0;
    int test_mass_count = 0;

    for (Long64_t i = 0; i < nEvents; i++) {
        ch->GetEntry(i);
        if (i % 10000 == 0) {
            cout << "Processing event (final correction) " << i << " / " << nEvents << "\r";
            cout.flush();
        }
        if (nclus < 2) continue;
        pass_clus_corr++;

        double best_dt = 1e9;
        int best_icl = -1, best_jcl = -1;
        best_pair_cuts(ecal_e, clus_nblk, ecal_x, ecal_y, clus_a_time, nclus, best_icl, best_jcl, best_dt);

        int imax1 = -1, imax2 = -1;
        double e1 = -1, e2 = -1;
        find_highest_energy(ecal_e, nclus, imax1, imax2, e1, e2);

        if (best_icl >= 0 && best_jcl >= 0 && ((best_icl == imax1 && best_jcl == imax2) || (best_icl == imax2 && best_jcl == imax1))) {
            pass_time_corr++;
            double e_corr[2] = {0., 0.}; // Stores corrected energy for each of the two clusters
            int idx_cl[2] = {best_icl, best_jcl};

            for (int ic = 0; ic < 2; ++ic) {
                int cl = idx_cl[ic];
                for (unsigned int b = 0; b < (unsigned int)ngoodADChits; ++b) {
                    int rawID = (int)goodblock_id[b];
                    int cid = (int)goodblock_cid[b];
                    auto it = blockID_to_idx.find(rawID);
                    // If blockID is known and belongs to the current cluster
                    if (it != blockID_to_idx.end() && cid == cl) {
                        e_corr[ic] += current_coeff[it->second] * goodblock_e[b];
                    }
                }
            }

            TVector3 pos1(ecal_x[best_icl], ecal_y[best_icl], z_calo);
            TVector3 pos2(ecal_x[best_jcl], ecal_y[best_jcl], z_calo);
            //TVector3 vertex(0, 0, z_target);
            TVector3 dir1 = (pos1 - vertex).Unit();
            TVector3 dir2 = (pos2 - vertex).Unit();
            TLorentzVector ph1_corr(dir1.X() * e_corr[0], dir1.Y() * e_corr[0], dir1.Z() * e_corr[0], e_corr[0]);
            TLorentzVector ph2_corr(dir2.X() * e_corr[1], dir2.Y() * e_corr[1], dir2.Z() * e_corr[1], e_corr[1]);
            double pi0_mass_corr = (ph1_corr + ph2_corr).M();
            Double_t opening_angle = dir1.Angle(dir2) * (180.0 / TMath::Pi());

            if (opening_angle < 3.5 || opening_angle > 8) continue;
            if (pi0_mass_corr > 0.02 && pi0_mass_corr < 0.4) {
                pass_mass_corr++;
                //h_pi0_mass_corr->Fill(pi0_mass_corr);
                test_mass_sum += pi0_mass_corr;
                test_mass_count++;
            }
        }
    }
    cout << "\n -----  Scaling coefficients to match PDG π⁰ mass.. -----" << endl;
    if (test_mass_sum > 0) {
        double mean_mass = test_mass_sum / test_mass_count;
        double scale_factor = pi0_mass_pdg / mean_mass;

        // Apply the scaling factor to the coefficients
        for (int i = 0; i < nblocks; ++i) {
            current_coeff[i] *= scale_factor;
        }
    }

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
                    e_corr[ic] += current_coeff[it->second] * goodblock_e[b];
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
            h_pi0_mass_corr->Fill(pi0_mass_corr);
        }
    }

    // -------------------- Plotting --------------------
    TCanvas *c = new TCanvas("c","Reconstructed #pi^{0} Invariant Mass Before and After Calibration", 800,600);
    h_pi0_mass->SetStats(0);
    h_pi0_mass_corr->SetStats(0);
    h_pi0_mass->GetXaxis()->SetTitle("M_{#pi^{0}} [GeV]");
    h_pi0_mass->GetYaxis()->SetTitle("Events");
    h_pi0_mass->SetLineWidth(2);
    h_pi0_mass_corr->SetLineWidth(2);
    h_pi0_mass->SetLineColor(kBlue);
    h_pi0_mass_corr->SetLineColor(kRed);
    h_pi0_mass_corr->Draw();
    h_pi0_mass->Draw("SAME");
    auto leg = new TLegend(0.6,0.7,0.9,0.9);
    leg->AddEntry(h_pi0_mass,      "Before calib", "l");
    //leg->AddEntry(h_pi0_mass_corr, "After calib",  "l");
    leg->AddEntry(h_pi0_mass_corr, ("After calib (Iter = " + to_string(n_iterations) + ")").c_str(), "l");
    leg->Draw();
    c->SaveAs("plots/ecal_pi0_mass_calib_iter.png");


    // Create a 2D histogram to visualize calibration coefficients in the detector geometry
    TH2F *h_coeff_map = new TH2F("h_coeff_map", "ECAL Block Calibration Coefficients;Column;Row",
        ncol, minCol, maxCol + 1l,
        nlin, minRow, maxRow + 1);

    // Fill the 2D histogram using the coeff[] array 
    for (int i = 0; i < nblocks; ++i) {
        int blockID = idx_to_blockID[i];
        int row = blockID_to_row[blockID];
        int col = blockID_to_col[blockID];
        double val = current_coeff[i];
        h_coeff_map->SetBinContent(col - minCol + 1, row - minRow + 1, val);  // ROOT bins start from 1
    }
    cout << "\n -----  Summary of the cuts -----" << endl;
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
    c_map->SaveAs("plots/ecal_coefficients_heatmap_iter.png");
}