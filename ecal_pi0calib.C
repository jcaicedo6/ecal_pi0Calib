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

void ecal_pi0calib() {
    const Double_t z_calo = 8.14; // position of calorimeter from the target in m
    const Double_t z_target = 0.09;    // position of target
    const Double_t z_origin = 0.0;
    const Double_t vertex_z = z_target - z_origin;  // position of vertex, where the pi0 is created, right in the middle of the target

    const int nbclusmax=100;// maximum number of clusters
    const int sizemax=1000;// maximum size of the cluster (how many blocks in each cluster)

    const double pi0_mass_pdg = 0.1349766;  // PDG pi0 mass in GeV
    TChain *ch = new TChain("T");

    string path = "/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign_lowcur_April28/rootfiles/";
    ifstream infile("lists_runs/runfiles_3173_to_3192.txt");
    string filename;

    while (getline(infile, filename)) {
        ch->Add((path + filename).c_str());
    }

    // Initialize the pointers before setting the branch address
    // Define the  energy, x, y of each cluster
    Double_t ecal_e[nbclusmax]; 
    Double_t ecal_x[nbclusmax];
    Double_t ecal_y[nbclusmax];

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
    Double_t clus_eblk[nbclusmax];       // block energy with highest energy in the cluster
    Double_t clus_idblk[nbclusmax];      // block id with highest energy in the cluster
    //Double_t energy_blk[1000];

    ch->SetBranchAddress("earm.ecal.clus.nblk", &clus_nblk);
    ch->SetBranchAddress("earm.ecal.clus.id", &clus_id);
    ch->SetBranchAddress("earm.ecal.clus.row", &clus_row);
    ch->SetBranchAddress("earm.ecal.clus.col", &clus_col);
    ch->SetBranchAddress("earm.ecal.clus.eblk", &clus_eblk);
    ch->SetBranchAddress("earm.ecal.idblk", &clus_idblk);
    //ch->SetBranchAddress("earm.ecal.a_p", &energy_blk);

    TH1F *h_pi0_mass = new TH1F("h_pi0_mass", "Pi0 Invariant Mass;M_{#pi^{0}} [GeV];Events", 80, 0, 0.6);
    TH1F *h_pi0_mass_corr = new TH1F("h_pi0_mass_reco", "Reconstructed Pi0 Invariant Mass;M_{#pi^{0}} [GeV];Events", 80, 0, 0.6);

    Long64_t nEvents = ch->GetEntries();
    cout << "Number of events: " << nEvents << endl;

    // --- Discover geometry from data: find min/max row & col ---
    int minRow = INT_MAX, maxRow = INT_MIN;
    int minCol = INT_MAX, maxCol = INT_MIN;

    for (Long64_t i = 0; i < nEvents; ++i) {
        ch->GetEntry(i);
        if (nclus < 1) continue;
        for (int cl = 0; cl < nclus && cl < nbclusmax; ++cl) {
            int r = static_cast<int>(clus_row[cl]);
            int c = static_cast<int>(clus_col[cl]);
            if (r < 0 || c < 0) continue;
            minRow = std::min(minRow, r);
            maxRow = std::max(maxRow, r);
            minCol = std::min(minCol, c);
            maxCol = std::max(maxCol, c);
        }
    }

    // Build geometry
    int nlin    = maxRow - minRow + 1;
    int ncol    = maxCol - minCol + 1;
    int nblocks = nlin * ncol;

    cout<<"Detected calo geometry: "<<ncol<<" cols × "<<nlin<<" rows"<<endl;
    cout << "Number of blocks: "<<nblocks<<endl;

    // ------------------- Determine the number of blocks -------------------

    // Mapping from original block ID to matrix index
    // Int_t num_btom[nblocks], num_mtob[nblocks];
    std::vector<int> num_btom(nblocks, -1), num_mtob(nblocks, -1);
    for(Int_t i=0;i<nblocks;i++){  

        num_btom[i]=-1;
        num_mtob[i]=-1;
    
    }

    int nbgood = 0;
    for (int i = 0; i < nblocks; ++i) {
        int ilin = i / ncol;
        int icol = i % ncol;
        if (icol != 0 && clus_eblk[0] > 0.0 && clus_eblk[1] > 0.0) {  // ignore bad column(s), e.g., 0
            num_mtob[nbgood] = i;
            num_btom[i] = nbgood;
            nbgood++;
        }
    }
    const int nblocksm = nbgood;

    cout << "Number of good blocks: " << nblocksm << endl;

    // Matrix for calibration
    TMatrixD A(nblocksm, nblocksm);
    //TMatrixD B(1, nblocksm);
    TVectorD B(nblocksm);
    A.Zero();
    B.Zero();
  
    // Block occupancy for later use
    std::vector<int> occupancy(nblocksm, 0);
    // indixes of good blocks
    std::vector<int> comp_idx(nblocks, -1);

    // Define the counters for the cuts
    Double_t pass_clus=0, pass_deltar=0, pass_time=0, pass_mass=0;
    // initialize the random number generator
    gRandom->SetSeed(0);
    //Double_t pi0_mass_pdg_random = pi0_mass_pdg;    

    for (Long64_t i = 0; i < nEvents; i++) {
        ch->GetEntry(i);
        
        // Check that there are at least two clusters
        if (nclus != 2) continue;
        pass_clus++;
        if ((ecal_e[0] + ecal_e[1]) < 0.5) continue;         // Minimum energy cut

        if (ecal_e[0] < 0.2 || ecal_e[1] < 0.2) continue;    // Minimum cluster energy cut Tune the 0.2 GeV threshold based on the noise level.
        if (clus_nblk[0] < 2 || clus_nblk[1] < 2) continue;  // Minimum number of blocks per cluster

        Double_t deltaR = sqrt(pow(ecal_x[0] - ecal_x[1], 2) + pow(ecal_y[0] - ecal_y[1], 2));
        if (deltaR < 0.07) continue;    // reject overlapping clusters
        pass_deltar++;

        if (clus_a_time[0] < 100 || clus_a_time[0] > 300) continue;  // Minimum time cut
        if (clus_a_time[1] < 100 || clus_a_time[1] > 300) continue;  // Minimum time cut
        if (fabs(clus_a_time[0] - clus_a_time[1]) > 10) continue;    // Maximum time difference between clusters
        pass_time++;
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

        // Calculate pi0 invariant mass
        Double_t pi0_mass = (ph1 + ph2).M();
        // Calculate opening angle between the two photons in degrees
        Double_t opening_angle = dir1.Angle(dir2) * (180.0 / TMath::Pi());
        // Apply cuts in the opening angle and pi0 mass
        if (opening_angle < 3.5 || opening_angle > 8) continue;  // The lower cut (e.g., < 6°) removes nearly collinear photon pairs → likely merged.
                                                                // The upper cut (e.g., > 80°) removes highly unphysical, possibly misreconstructed pairs.
        if (pi0_mass <= 0.0 || pi0_mass >= 0.4) continue;    // Making a cut on the pi0 mass, e.g., between 0.06 and 0.6 GeV

        
        // Determine the scale factor
        // avoid zero division
        if (pi0_mass <= 0) continue;
        pass_mass++;

        // Filling the uncorrected pi0 mass histogram
        h_pi0_mass->Fill(pi0_mass);

        // compute expected total π0 energy as using the scale factor
        //double pi0_mass_smeared = gRandom->Gaus(pi0_mass_pdg, 0.001); // 5 MeV width max
        double expected_E = (ecal_e[0] + ecal_e[1]) * (pi0_mass_pdg / pi0_mass);
        //if (expected_E / (ecal_e[0] + ecal_e[1]) > 3.0) continue;

        
        // **Loop over all blocks in cluster **:
        // zero the per-block energy accumulator
        std::vector<double> energy(nblocksm, 0.0);
    
        int idx_blk = 0;
  
        // Loop over all clusters and blocks per event
        for (int cl = 0; cl < nclus; ++cl) {
            // Energy of block with highest energy in cluster
            double e_main_blk = clus_eblk[cl];
            // Number of blocks in cluster
            int nblk = clus_nblk[cl];
            // remaining energy to be shared among other blocks
            double e_rem = (nblk > 1) ? (ecal_e[cl] - e_main_blk) / (nblk - 1) : 0.0;
            // Loop over all blocks in cluster
            for (int b = 0; b < nblk; ++b, ++idx_blk) {
                //int rawID = row*ncol + col;
                // scaning rawID in clus_id
                int rawID = static_cast<int>(clus_id[idx_blk]);
                // skip bad blocks
                if (rawID<0 || rawID>=nblocks) continue;
                // Matching rawID to num_btom to find block id
                int iblock   = num_btom[rawID];
                // save block id
                comp_idx[idx_blk] = iblock;
                if (iblock<0 || iblock>=nblocksm) continue;  // skip bad blocks

                // Matching rawID to clus_idblk to find max energy block id and assigning there the max energy
                if (rawID == (int)clus_idblk[cl]) {
                    energy[iblock] += e_main_blk;  
                }
                else {
                    energy[iblock] += e_rem;
                }   
                if (energy[iblock] > 0.0) {
                    occupancy[iblock]++;      // update occupancy
                }
                // Commented section for debugging
                //std::cout << "Cluster " << cl << ", Block " << b
                  //  << ", rawID: " << rawID << ", iblock: " << iblock
                  //  << ", energy: " << energy[iblock] << std::endl;

                
            }    
            //std::cout << "[DEBUG] clus_nblk = " << clus_nblk[cl]
              //      << " ecal_e = " << ecal_e[cl]
              //      << " clus_eblk = " << clus_eblk[cl]
              //      << " computed e_rem = " << e_rem << std::endl;        
            
        }

        //int active = 0;
        //for (double e : energy) if (e!=0.) ++active;
          //  if (active==0) {
          //  std::cout<<"[DEBUG] Event "<<i<<" passed cuts but energy[] all zero\n";

        //}
        
        // --------------- Debugging --------------------
        //print the sum of energy in each block in cluster and comparing it with the ecal_e
        /*Double_t energy_sum_clus1 = 0, energy_sum_clus2 = 0;
        for (int b = 0; b < clus_nblk[0]; ++b) {
            int row   = clus_row[b];
            int col   = clus_col[b];
            //int rawID = row*ncol + col;
            int rawID = static_cast<int>(clus_id[b]);
            if (rawID<0 || rawID>=nblocks) continue;
            int idx   = num_btom[rawID];
            if (idx<0 || idx>=nblocksm) continue;
            energy_sum_clus1 += energy[idx];
        }
        for (int b = 0; b < clus_nblk[1]; ++b) {
            int row   = clus_row[b];
            int col   = clus_col[b];
            //int rawID = row*ncol + col;
            int rawID = static_cast<int>(clus_id[b]);
            if (rawID<0 || rawID>=nblocks) continue;
            int idx   = num_btom[rawID];
            if (idx<0 || idx>=nblocksm) continue;
            energy_sum_clus2 += energy[idx];
        }*/

        
        //std::cout<<"[DEBUG] Event "<<i<<" energy_sum_clus1 = "<<energy_sum_clus1<<" ecal_e[0] = "<<ecal_e[0]<<" energy_sum_clus2 = "<<energy_sum_clus2<<" ecal_e[1] = "<<ecal_e[1]<<"\n";

        // Now fill A and B:
        for (int j = 0; j < nblocksm; ++j) {
            for (int k = 0; k < nblocksm; ++k) {
                A(j,k) += energy[j] * energy[k];
            }
            //B(0,j) += energy[j] * expected_E;
            B(j) += energy[j] * expected_E;
        }           
    }
    cout<<"Events passing number of clusters cut: "<<pass_clus<<endl;
    cout<<"Events passing deltaR cut: "<<pass_deltar<<endl;
    cout<<"Events passing time cut: "<<pass_time<<endl;
    cout<<"Events passing all cuts: "<<pass_mass<<endl;
  


    // ===================== OCCUPANCY PRUNING & REGULARIZATION =====================
    const int N_min = 1;
    for (int i = 0; i < nblocksm; ++i) {
        if (occupancy[i] < N_min) {
            // zero out row and column i
            for (int j = 0; j < nblocksm; ++j) {
                A(i,j) = A(j,i) = 0.0;
            }
            //B(0,i) = 0.0;
            B(i) = 0.0;
        }
    }
    // add small diagonal term λ·A_ii
    const double lambda = 1e-2;
    for (int i = 0; i < nblocksm; ++i) {
        A(i,i) += lambda * A(i,i);
    }
        
    // -------------------- Invert the Matrix --------------------
    //A.Print();
    std::cout
        <<"[DEBUG] A(0,0)="<<A(0,0)
        <<", B(0,0)="<<B(0)
        <<"\n";
        for (int i = 0; i < std::min(5,nblocksm); ++i) {
        std::cout<<"  B(0,"<<i<<")="<<B(i)<<"\n";
        }

    //Double_t coeff[nblocks];
    std::vector<double> coeff(nblocks, 1.0);


    //for(Int_t i=0;i<nblocks;i++){coeff[i]=1.;} // Set all coefficients to 1 for missing blocks
    // overwrite good blocks

    // Check if the matrix is invertible
    // invert & solve C=B·A⁻¹
    /*if (TMath::Abs(A.Determinant()) < 1e-10) {
        cout << "Matrix is singular. Cannot invert." << endl;
     
    } else {
        A.Invert();
        // Continue with the inversion and other calculations
    }*/

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
    for (int compID = 0; compID < nblocksm; ++compID) {
        int rawID = num_mtob[compID];     // reverse map
        //double c = C(0,compID);
        double c = C(compID);

        if (occupancy[compID] < N_min || c < 1e-2 || c > 10) {
            coeff[rawID] = 1.0;
        } else {
            coeff[rawID] = c;
        //coeff[rawID] = C(0, compID);      // use same compID as above
        }
    }
        // now overwrite just the “good” ones:
        //for (int i = 0; i < nblocksm; ++i) {
        //  int rawID = num_mtob[i];  // maps 0…nblocksm-1 back to 0…nblocks-1
            //double c = C(0,i);
        //  coeff[ rawID ] = C(0, i);

            //if (occupancy[i] < N_min || c < 1e-2 || c > 10) {
            //  coeff[rawID] = 1.0;
            //} else {
            //   coeff[rawID] = c;
            //}
        //}
    

    // -------------------- Saving the coefficients in a text file --------------------
    cout<<"Writing coefficients to file..."<<endl;

    ofstream coeff_file("ecal_block_calibration_factors.txt");
    coeff_file << "# BlockID\tCalibrationCoeff\n";

    for (int i = 0; i < nblocks; i++) {
        coeff_file << i << "\t" << coeff[i] << endl;
        cout << "BlockID: " << i << "\tCalibrationCoeff: " << coeff[i] << endl;
    }

    coeff_file.close();

    // Now, with the coefficients in hand, I can use them to calibrate the Ecal energy in the code.
    // Use the coefficients to calibrate the Ecal energy and plot the pi0 invariant mass
    // with and without the correction

    // counters for passing events
    Double_t pass_clus_corr=0, pass_deltar_corr=0, pass_time_corr=0, pass_mass_corr=0;
 
    for (Long64_t i = 0; i < nEvents; i++) {
        ch->GetEntry(i);
        if (nclus != 2) continue;
        if ((ecal_e[0] + ecal_e[1]) < 0.5) continue;    // reject low energy events
        pass_clus_corr++;

        if (ecal_e[0] < 0.2 || ecal_e[1] < 0.2) continue;    // Minimum cluster energy cut Tune the 0.2 GeV threshold based on your noise level.
        if (clus_nblk[0] < 2 || clus_nblk[1] < 2) continue;  // Minimum number of blocks per cluster

        Double_t deltaR = sqrt(pow(ecal_x[0] - ecal_x[1], 2) + pow(ecal_y[0] - ecal_y[1], 2));
        if (deltaR < 0.07) continue; // reject overlapping clusters
        pass_deltar_corr++;

        if (clus_a_time[0] < 100 || clus_a_time[0] > 300) continue;  // reject events with bad timing
        if (clus_a_time[1] < 100 || clus_a_time[1] > 300) continue;  // reject events with bad timing
        if (fabs(clus_a_time[0] - clus_a_time[1]) > 10) continue;   // reject events with bad timing
        pass_time_corr++;

        //cout << "Processing event " << i << "\r";
        if (i % 1000 == 0) cout << "Processing event " << i << " / " << nEvents << endl;



        // **After** correction: rebuild each photon’s energy by summing
        //     per‑block energies * coefficients
        double e_corr[2] = {0., 0.};
        int id_blk = 0;

        for (int cl =0; cl < nclus; ++cl) {
            // Energy of block with highest energy in cluster
            double e_main_blk = clus_eblk[cl];
            // Number of blocks in cluster
            int nblk = clus_nblk[cl];
            // remaining energy to be shared among other blocks
            double e_rem = (nblk > 1) ? (ecal_e[cl] - e_main_blk) / (nblk - 1) : 0.0;

            // Loop over all blocks in cluster
            for (int b = 0; b < nblk; ++b, ++id_blk) {
                int rawID = static_cast<int>(clus_id[id_blk]);
                if (rawID<0 || rawID>=nblocks) continue;
                int iblock   = num_btom[rawID];
                if (iblock<0 || iblock>=nblocksm) continue;

                if (rawID == (int)clus_idblk[cl]) {
                    e_corr[cl] += coeff[rawID] * e_main_blk;
                } else {
                    e_corr[cl] += coeff[rawID] * e_rem;
                }
            }

        }

        TVector3 pos1(ecal_x[0], ecal_y[0], z_calo);  // in m
        TVector3 pos2(ecal_x[1], ecal_y[1], z_calo);  // in m

        // uniform smearing ±15 cm around (0,0,0.09)
        double vertex_x_smeared = gRandom->Uniform(-0.15, +0.15);
        double vertex_y_smeared = gRandom->Uniform(-0.15, +0.15);
        double vertex_z_smeared = z_target + gRandom->Uniform(-0.15, +0.15);
        TVector3 vertex(vertex_x_smeared, vertex_y_smeared, vertex_z_smeared);
        TVector3 dir1 = (pos1 - vertex).Unit();
        TVector3 dir2 = (pos2 - vertex).Unit();

        // rebuild corrected TLorentzVectors
        TLorentzVector ph1_corr(dir1.X() * e_corr[0], dir1.Y() * e_corr[0], dir1.Z() * e_corr[0], e_corr[0]);
        TLorentzVector ph2_corr(dir2.X() * e_corr[1], dir2.Y() * e_corr[1], dir2.Z() * e_corr[1], e_corr[1]);
        
        double pi0_mass_corr = (ph1_corr + ph2_corr).M();

        Double_t opening_angle = dir1.Angle(dir2) * (180.0 / TMath::Pi());
        if (opening_angle < 3.5 || opening_angle > 8) continue;  // The lower cut (e.g., < 6°) removes nearly collinear photon pairs → likely merged.
                                                                // The upper cut (e.g., > 80°) removes highly unphysical, possibly misreconstructed pairs.

        if (pi0_mass_corr <= 0.0 || pi0_mass_corr >= 0.4) continue;    // Making a cut on the pi0 mass, e.g., between 0.06 and 0.6 GeV

        // fill histograms
        pass_mass_corr++;
        h_pi0_mass_corr->Fill(pi0_mass_corr);
        
       
    }
    double mu = h_pi0_mass_corr->GetMean(); 
    double s = pi0_mass_pdg / mu;
    for (int raw = 0; raw < nblocks; ++raw)
        coeff[raw] *= s;

    h_pi0_mass_corr->Reset();  // clear previous content

    for (Long64_t i = 0; i < nEvents; i++) {
        ch->GetEntry(i);
        if (nclus != 2) continue;
        if ((ecal_e[0] + ecal_e[1]) < 0.5) continue;    // reject low energy events

        if (ecal_e[0] < 0.2 || ecal_e[1] < 0.2) continue;    // Minimum cluster energy cut Tune the 0.2 GeV threshold based on your noise level.
        if (clus_nblk[0] < 2 || clus_nblk[1] < 2) continue;  // Minimum number of blocks per cluster

        Double_t deltaR = sqrt(pow(ecal_x[0] - ecal_x[1], 2) + pow(ecal_y[0] - ecal_y[1], 2));
        if (deltaR < 0.07) continue; // reject overlapping clusters

        if (clus_a_time[0] < 100 || clus_a_time[0] > 300) continue;  // reject events with bad timing
        if (clus_a_time[1] < 100 || clus_a_time[1] > 300) continue;  // reject events with bad timing
        if (fabs(clus_a_time[0] - clus_a_time[1]) > 10) continue;   // reject events with bad timing

        //cout << "Processing event " << i << "\r";
        if (i % 1000 == 0) cout << "Processing event " << i << " / " << nEvents << endl;

        // **After** correction: rebuild each photon’s energy by summing
        //     per‑block energies * coefficients
        double e_corr[2] = {0., 0.};
        int id_blk = 0;

        for (int cl =0; cl < nclus; ++cl) {
            // Energy of block with highest energy in cluster
            double e_main_blk = clus_eblk[cl];
            // Number of blocks in cluster
            int nblk = clus_nblk[cl];
            // remaining energy to be shared among other blocks
            double e_rem = (nblk > 1) ? (ecal_e[cl] - e_main_blk) / (nblk - 1) : 0.0;

            // Loop over all blocks in cluster
            for (int b = 0; b < nblk; ++b, ++id_blk) {
                int rawID = static_cast<int>(clus_id[id_blk]);
                if (rawID<0 || rawID>=nblocks) continue;
                int iblock   = num_btom[rawID];
                if (iblock<0 || iblock>=nblocksm) continue;

                if (rawID == (int)clus_idblk[cl]) {
                    e_corr[cl] += coeff[rawID] * e_main_blk;
                } else {
                    e_corr[cl] += coeff[rawID] * e_rem;
                }
            }

        }

        TVector3 pos1(ecal_x[0], ecal_y[0], z_calo);  // in m
        TVector3 pos2(ecal_x[1], ecal_y[1], z_calo);  // in m

        // uniform smearing ±15 cm around (0,0,0.09)
        double vertex_x_smeared = gRandom->Uniform(-0.15, +0.15);
        double vertex_y_smeared = gRandom->Uniform(-0.15, +0.15);
        double vertex_z_smeared = z_target + gRandom->Uniform(-0.15, +0.15);
        TVector3 vertex(vertex_x_smeared, vertex_y_smeared, vertex_z_smeared);
        TVector3 dir1 = (pos1 - vertex).Unit();
        TVector3 dir2 = (pos2 - vertex).Unit();

        // rebuild corrected TLorentzVectors
        TLorentzVector ph1_corr(dir1.X() * e_corr[0], dir1.Y() * e_corr[0], dir1.Z() * e_corr[0], e_corr[0]);
        TLorentzVector ph2_corr(dir2.X() * e_corr[1], dir2.Y() * e_corr[1], dir2.Z() * e_corr[1], e_corr[1]);
        
        double pi0_mass_corr = (ph1_corr + ph2_corr).M();

        Double_t opening_angle = dir1.Angle(dir2) * (180.0 / TMath::Pi());
        if (opening_angle < 3.5 || opening_angle > 8) continue;  // The lower cut (e.g., < 6°) removes nearly collinear photon pairs → likely merged.
                                                                // The upper cut (e.g., > 80°) removes highly unphysical, possibly misreconstructed pairs.

        if (pi0_mass_corr <= 0.0 || pi0_mass_corr >= 0.4) continue;    // Making a cut on the pi0 mass, e.g., between 0.06 and 0.6 GeV

        // fill histograms
        h_pi0_mass_corr->Fill(pi0_mass_corr);
        
       
    }
    // -------------------- Plotting --------------------
    TCanvas *c = new TCanvas("c","before/after", 800,600);
    h_pi0_mass->SetLineColor(kRed);
    h_pi0_mass_corr->SetLineColor(kBlue);
    h_pi0_mass->Draw();
    h_pi0_mass_corr->Draw("SAME");
    auto leg = new TLegend(0.6,0.7,0.9,0.9);
    leg->AddEntry(h_pi0_mass,      "Before calib", "l");
    leg->AddEntry(h_pi0_mass_corr, "After calib",  "l");
    leg->Draw();


    // Create a 2D histogram to visualize calibration coefficients
    TH2F *h_coeff_map = new TH2F("h_coeff_map", "ECAL Block Calibration Coefficients;Column;Row",
        ncol, minCol, minCol + ncol,
        nlin, minRow, minRow + nlin);

    // Fill the 2D histogram using the coeff[] array
    for (int row = 0; row < nlin; ++row) {
        for (int col = 0; col < ncol; ++col) {
            int rawID = row * ncol + col;
            if (rawID < 0 || rawID >= nblocks) continue;

            double val = coeff[rawID];
            h_coeff_map->SetBinContent(col + 1, row + 1, val);  // ROOT bins start from 1
        }
    }
    cout<<"Events passing number of clusters cut after calib: "<<pass_clus_corr<<endl;
    cout<<"Events passing deltaR cut after calib: "<<pass_deltar_corr<<endl;
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

    //TCanvas *c_pi0InvM = new TCanvas("c_pi0InvM", "c_pi0InvM", 800, 600);
    //h_pi0_mass->Draw();

    //TCanvas *c_photons_time = new TCanvas("c_photons_time", "c_photons_time", 800, 600);
    //h_photons_time->Draw("colz");



    //TFile *f = new TFile("pi0_mass.root", "RECREATE");
    //h_pi0_mass->Write();
    //f->Close();
}