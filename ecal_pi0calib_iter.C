#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TDecompSVD.h"
#include "TRandom3.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

// This script calibrates the energy of the ecal clusters in a iterative way
// run with “root -l ecal_pi0calib_iter.C++”
void ecal_pi0calib_iter() {
    // constants
    const Double_t z_calo    = 8.14;                         // position of calorimeter from the target in m
    const Double_t z_target  = 0.09;                         // position of target
    const Double_t pi0_mass_pdg = 0.1349766;                 // PDG pi0 mass in GeV
    const int    nbclusmax  = 100;                           // maximum number of clusters

    // Build chain
    TChain *ch = new TChain("T");
    ifstream infile("lists_runs/runfiles_3173_to_3192.txt");
    string filename;
    while (getline(infile, filename))
        ch->Add((string("/volatile/halla/sbs/sbs-gep/GEP_REPLAYS/GEP1/LH2/prod_realign_lowcur_April28/rootfiles/")+filename).c_str());

    // Branches
    Double_t ecal_e[nbclusmax], ecal_x[nbclusmax], ecal_y[nbclusmax];
    Double_t nclus, clus_a_time[1000];
    Double_t clus_nblk[nbclusmax], clus_id[nbclusmax],
             clus_row[nbclusmax], clus_col[nbclusmax],
             clus_eblk[nbclusmax], clus_idblk[nbclusmax];

    ch->SetBranchAddress("earm.ecal.clus.e",      &ecal_e);            // ecal cluster energy
    ch->SetBranchAddress("earm.ecal.clus.x",      &ecal_x);            // ecal cluster x position
    ch->SetBranchAddress("earm.ecal.clus.y",      &ecal_y);            // ecal cluster y position
    ch->SetBranchAddress("earm.ecal.nclus",       &nclus);             // number of clusters
    ch->SetBranchAddress("earm.ecal.clus.atimeblk",&clus_a_time);      // cluster arrival time
    ch->SetBranchAddress("earm.ecal.clus.nblk",   &clus_nblk);         // number of blocks
    ch->SetBranchAddress("earm.ecal.clus.id",     &clus_id);           // cluster id
    ch->SetBranchAddress("earm.ecal.clus.row",    &clus_row);          // cluster row
    ch->SetBranchAddress("earm.ecal.clus.col",    &clus_col);          // cluster column
    ch->SetBranchAddress("earm.ecal.clus.eblk",   &clus_eblk);         // energy of each block
    ch->SetBranchAddress("earm.ecal.idblk",       &clus_idblk);        // id of each block

    Long64_t nEvents = ch->GetEntries();
    cout << "Events: " << nEvents << "\n";

    // --- Discover geometry from data: find min/max row & col ---
    int minRow=INT_MAX, maxRow=INT_MIN,
        minCol=INT_MAX, maxCol=INT_MIN;
    for (Long64_t i=0; i<nEvents; ++i) {
        ch->GetEntry(i);
        if (nclus<1) continue;
        for (int cl=0; cl<nclus && cl<nbclusmax; ++cl) {
            int r = clus_row[cl], c = clus_col[cl];
            if (r<0||c<0) continue;  // ignore bad column(s), e.g., 0
            minRow = min(minRow,r);
            maxRow = max(maxRow,r);
            minCol = min(minCol,c);
            maxCol = max(maxCol,c);
        }
    }
    int nlin = maxRow-minRow+1,
        ncol = maxCol-minCol+1,
        nblocks = nlin*ncol;
    
    cout<<"Detected calo geometry: "<<ncol<<" cols × "<<nlin<<" rows"<<endl;
    cout << "Number of blocks: "<<nblocks<<endl;

    // ------------------- Determine the number of blocks -------------------
    // Mapping from original block ID to matrix index AND vice versa
    vector<int> num_btom(nblocks,-1), num_mtob(nblocks,-1);
    int nbgood=0;
    for (int i=0;i<nblocks;++i) {
        int col = i % ncol;
        if (col!=0 && clus_eblk[0] > 0.0 && clus_eblk[1] > 0.0) { // ignore bad column(s), e.g., 0
            num_btom[i]=nbgood;
            num_mtob[nbgood]=i;
            ++nbgood;
        }
    }
    const int nblocksm = nbgood;

    cout << "Number of good blocks: " << nblocksm << endl;

    // Initialize coefficients = 1 for now
    vector<double> coeff(nblocks,1.0);
    vector<double> old_coeff = coeff;        // make a copy

    // Histograms for final comparison
    TH1F *h_before = new TH1F("h_before",
      "Pi^{0} Mass Before/After;M_{#pi^{0}}[GeV];Events",90,0,0.6);
    TH1F *h_after = (TH1F*)h_before->Clone("h_after");
    h_after->SetTitle("Corrected");

    gRandom->SetSeed(0);

    // Iterative calibration
    const int nIter = 100;          // number of passes
    for (int iter=0; iter<nIter; ++iter) {
        cout<<"Starting the "<<iter+1<<" Iteration"<<endl;

        // copy old coefficients
        old_coeff = coeff;
        // zero matrices
        TMatrixD A(nblocksm,nblocksm); A.Zero();
        TVectorD B(nblocksm);          B.Zero();
        // Block occupancy
        vector<int> occupancy(nblocksm,0);

        // 
        for (Long64_t i=0; i<nEvents; ++i) {
            ch->GetEntry(i);
            // Check if there are two clusters
            if (nclus!=2) continue;
            if (ecal_e[0]+ecal_e[1]<0.5) continue;               // Minimum cluster energy
            if (ecal_e[0] < 0.2 || ecal_e[1] < 0.2) continue;    // Minimum cluster energy cut Tune the 0.2 GeV threshold based on the noise level.
            if (clus_nblk[0] < 2 || clus_nblk[1] < 2) continue;  // Minimum number of blocks per cluster
            // Cuts on cluster energy, size, timing, etc.
            double dR = sqrt(pow(ecal_x[0]-ecal_x[1],2)+
                             pow(ecal_y[0]-ecal_y[1],2));
            if (dR<0.07) continue;
            if (clus_a_time[0] < 100 || clus_a_time[0] > 300) continue;  // Minimum time cut
            if (clus_a_time[1] < 100 || clus_a_time[1] > 300) continue;  // Minimum time cut
            if (fabs(clus_a_time[0]-clus_a_time[1])>10) continue;        // Maximum time difference between clusters

            if (i % 500 == 0) cout << "Processing event " << i << " / " << nEvents << " for the " << iter+1 << " Iteration" <<  endl;

            // Compute pi0 mass and expected total E
            // Uniform smearing of the vertex within ±15 cm around (0,0,0.09)
            double vx_smeared = gRandom->Uniform(-0.15, +0.15);
            double vy_smeared = gRandom->Uniform(-0.15, +0.15);
            double vz_smeared = z_target + gRandom->Uniform(-0.15, +0.15);
            // compute the vertex vectors
            TVector3 v1(ecal_x[0],ecal_y[0],z_calo),
                     v2(ecal_x[1],ecal_y[1],z_calo),
                     vz(vx_smeared,vy_smeared,vz_smeared);
            // compute the direction vectors
            TVector3 d1=(v1-vz).Unit(), d2=(v2-vz).Unit();
            // compute the Lorentz vectors of the clusters
            TLorentzVector P1(d1*ecal_e[0],ecal_e[0]),
                            P2(d2*ecal_e[1],ecal_e[1]);

            // Compute expected opening angle between the photons in degrees
            double op_angle = d1.Angle(d2) * (180.0 / TMath::Pi());
            // Cuts on opening angle 
            if (op_angle < 3.5 || op_angle > 8) continue;  // The lower cut (e.g., < 3.5°) removes nearly collinear photon pairs → likely merged.
                                                           // The upper cut (e.g., > 8°) removes highly unphysical, possibly misreconstructed pairs.
            // Compute expected pi0 mass
            double m = (P1+P2).M();
            if (m <= 0.0 || m >= 0.4) continue;
            double Eexp = (ecal_e[0]+ecal_e[1])*(pi0_mass_pdg/m);

            // **Loop over all blocks in cluster **:
            // zero the per-block energy accumulator
            vector<double> Eblock(nblocksm,0.);
            int idx=0;
            // Loop over all clusters and blocks per event
            for (int cl=0;cl<2;++cl) { 
                double emax = clus_eblk[cl];                                        // Energy of block with highest energy in cluster
                int    nbl = clus_nblk[cl];                                         // Number of blocks in cluster
                double er = (nbl>1?(ecal_e[cl]-emax)/(nbl-1):0.);                   // remaining energy to be shared among other blocks
                // Loop over all blocks in cluster
                for (int b=0;b<nbl;++b,++idx) {
                    int raw = clus_id[idx];                                         // scaning rawID in clus_id
                    if (raw<0||raw>=nblocks) continue;                              // check if rawID is valid
                    int ib = num_btom[raw];                                         // convert rawID to blockID, Matching rawID to num_btom to find block id
                    if (ib<0||ib>=nblocksm) continue;
                    double share = (raw == (int)clus_idblk[cl] ? emax : er);        // find the block id with highest energy and assigning there the max energy
                    Eblock[ib] += share * coeff[ raw ];                             // add the energy to the block and scale by coeff
                    occupancy[ib]++;                                                // increment occupancy
                }
            }

            // Now fill A and B:
            for (int j=0;j<nblocksm;++j) {
                for (int k=0;k<nblocksm;++k)
                    A(j,k) += Eblock[j]*Eblock[k];
                B(j) += Eblock[j]*Eexp;
            }
        }

        // ===================== OCCUPANCY PRUNING & REGULARIZATION =====================
        int N_ocmin = 1;                                   // minimum occupancy
        const double lambda = 1e-2;                        // regularization
        const double α      = 0.1;                         // relaxing update factor for coeff, it smooths convergence. Lower α (e.g. to 0.3) to further damp oscillations.
        for (int i=0;i<nblocksm;++i) {
            if (occupancy[i]<N_ocmin) {
                for (int j=0;j<nblocksm;++j)
                    A(i,j)=A(j,i)=0;
                B(i)=0;
            }
            A(i,i) += lambda*A(i,i);                       // regularization, add small diagonal term λ·A_ii
        }

        // SVD decomposition of A and solve C = A⁺·B (see https://root.cern.ch/doc/master/classTDecompSVD.html)
        // the coefficients of the linear system 
        // C = A^-1 * B where A^-1 is the inverse of A, B is the vector of unknowns, and C is the vector of calibration coefficients. originally A * C = B
        TDecompSVD svd(A);
        if (!svd.Decompose()) {
            cerr<<"SVD failed at iter "<<iter<<"\n";
            break;
        }
        // Prepare solution vector
        // Copy B into C, then solve C = A⁺·B
        TVectorD C = B;             // initialize C with RHS
        if (!svd.Solve(C)) {        // Solve in place: C ← pseudo-inverse(A)·C
            cerr<<"Solve failed at iter "<<iter<<"\n";
            break;
        }

        // update coeff[]
        /*for (int ib=0;ib<nblocksm;++ib) {
            double c = C(ib);
            int raw = num_mtob[ib];
            if (occupancy[ib]<N_ocmin || c<0.01 || c>10) {
                coeff[raw] = 1.0;
            } else {
                coeff[raw] = c;
            }
        }*/
        double delta2 = 0;
        for (int ib=0; ib<nblocksm; ++ib) {
            double newc = C(ib);
            int raw     = num_mtob[ib];
            if (occupancy[ib] < N_ocmin || newc < 0.01 || newc > 10) {
              newc = 1.0;
            } else{
                if (iter > 0) {
                    coeff[raw] = coeff[raw] * (1.0 - α) + newc * α;
                } else {
                    coeff[raw] = newc;
                }
            }
            // compute the RMS coefficients shift
            double d = coeff[raw] - old_coeff[raw];
            delta2 += d * d;
        }
        // ...existing iterative calibration code...
        cout << "Iteration " << iter + 1 << " shift RMS = " << sqrt(delta2 / nblocksm) << "\n";
        if (sqrt(delta2 / nblocksm) < 5e-2) {
            cout << "Converged after " << iter + 1 << " iterations.\n";
            break;
        }
    }

    // Compute the corrected mass (to calculate the scaling factor)
    double sum_corrected_mass = 0.0;
    int count_corrected_mass = 0;

    for (Long64_t i=0;i<nEvents;++i) {
        ch->GetEntry(i);
        if (nclus!=2) continue;
        if (ecal_e[0]+ecal_e[1]<0.5) continue;
        if (ecal_e[0]<0.2||ecal_e[1]<0.2) continue;
        if (clus_nblk[0] < 2 || clus_nblk[1] < 2) continue;  // Minimum number of blocks per cluster
        double dR = sqrt(pow(ecal_x[0]-ecal_x[1],2)+
                         pow(ecal_y[0]-ecal_y[1],2));
        if (dR<0.07) continue;
        if (clus_a_time[0] < 100 || clus_a_time[0] > 300) continue;  // Minimum time cut
        if (clus_a_time[1] < 100 || clus_a_time[1] > 300) continue;  // Minimum time cut
        if (fabs(clus_a_time[0]-clus_a_time[1])>10) continue;

        if (i % 500 == 0) cout << "Processing event " << i << " / " << nEvents << " for applying the calibration" <<  endl;

        // raw mass
        // uniform smearing ±15 cm around (0,0,0.09)
        double vx_smeared = gRandom->Uniform(-0.15, +0.15);
        double vy_smeared = gRandom->Uniform(-0.15, +0.15);
        double vz_smeared = z_target + gRandom->Uniform(-0.15, +0.15);
        TVector3 v1(ecal_x[0],ecal_y[0],z_calo),
                 v2(ecal_x[1],ecal_y[1],z_calo),
                 vz(vx_smeared,vy_smeared,vz_smeared);
        TVector3 d1=(v1-vz).Unit(), d2=(v2-vz).Unit();
        TLorentzVector P1(d1*ecal_e[0],ecal_e[0]),
                        P2(d2*ecal_e[1],ecal_e[1]);

        double op_ang = d1.Angle(d2) * (180.0 / TMath::Pi());
        if (op_ang<3.5||op_ang>8) continue;
        double m_raw = (P1+P2).M();
        if (m_raw>=0.0&&m_raw<=0.4) h_before->Fill(m_raw);

        // corrected energy
        double Ecor[2]={0,0}; int idx=0;
        for (int cl=0;cl<2;++cl) {
            double emax = clus_eblk[cl];
            int    nbl = clus_nblk[cl];
            double er = (nbl>1?(ecal_e[cl]-emax)/(nbl-1):0.);
            for (int b=0;b<nbl;++b,++idx) {
                int raw=clus_id[idx];
                if (raw<0||raw>=nblocks) continue;
                double share=(raw==(int)clus_idblk[cl]?emax:er);
                Ecor[cl]+=share*coeff[raw];
            }
        }
        TLorentzVector Q1(d1*Ecor[0],Ecor[0]),
                        Q2(d2*Ecor[1],Ecor[1]);
        double m_cor = (Q1+Q2).M();
        if (m_cor >= 0.0 && m_cor <= 0.4) {
            sum_corrected_mass += m_cor;
            count_corrected_mass++;
        }
    }

    // Calculate the scaling factor
    cout << "Scaling coefficients to match PDG π⁰ mass..." << endl;
    if (count_corrected_mass > 0) {
        double mean_corrected_mass = sum_corrected_mass / count_corrected_mass;
        double scale_factor = pi0_mass_pdg / mean_corrected_mass;

        // Apply the scaling factor to the coefficients
        for (int raw = 0; raw < nblocks; ++raw) {
            coeff[raw] *= scale_factor;
        }
    }
    // Final histogram filling with scaled coefficients
    for (Long64_t i = 0; i < nEvents; ++i) {
        ch->GetEntry(i);
        if (nclus != 2) continue;
        if (ecal_e[0] + ecal_e[1] < 0.5) continue;
        if (ecal_e[0] < 0.2 || ecal_e[1] < 0.2) continue;
        if (clus_nblk[0] < 2 || clus_nblk[1] < 2) continue;
    
        double dR = sqrt(pow(ecal_x[0] - ecal_x[1], 2) + pow(ecal_y[0] - ecal_y[1], 2));
        if (dR < 0.07) continue;
        if (clus_a_time[0] < 100 || clus_a_time[0] > 300) continue;
        if (clus_a_time[1] < 100 || clus_a_time[1] > 300) continue;
        if (fabs(clus_a_time[0] - clus_a_time[1]) > 10) continue;
    
        double vx_smeared = gRandom->Uniform(-0.15, +0.15);
        double vy_smeared = gRandom->Uniform(-0.15, +0.15);
        double vz_smeared = z_target + gRandom->Uniform(-0.15, +0.15);
    
        TVector3 v1(ecal_x[0], ecal_y[0], z_calo), v2(ecal_x[1], ecal_y[1], z_calo),
                 vz(vx_smeared, vy_smeared, vz_smeared);
    
        TVector3 d1 = (v1 - vz).Unit(), d2 = (v2 - vz).Unit();
        TLorentzVector P1(d1 * ecal_e[0], ecal_e[0]), P2(d2 * ecal_e[1], ecal_e[1]);
    
        double op_angle = d1.Angle(d2) * (180.0 / TMath::Pi());
        if (op_angle < 3.5 || op_angle > 8) continue;
    
        double Ecor[2] = {0, 0};
        int idx = 0;
        for (int cl = 0; cl < 2; ++cl) {
            double emax = clus_eblk[cl];
            int nbl = clus_nblk[cl];
            double er = (nbl > 1 ? (ecal_e[cl] - emax) / (nbl - 1) : 0.);
            for (int b = 0; b < nbl; ++b, ++idx) {
                int raw = clus_id[idx];
                if (raw < 0 || raw >= nblocks) continue;
                double share = (raw == (int)clus_idblk[cl] ? emax : er);
                Ecor[cl] += share * coeff[raw];
            }
        }
    
        TLorentzVector Q1(d1 * Ecor[0], Ecor[0]), Q2(d2 * Ecor[1], Ecor[1]);
        double m_cor = (Q1 + Q2).M();
        if (m_cor >= 0.0 && m_cor <= 0.4) h_after->Fill(m_cor);
    }

    cout<<"-------------------- Difference between before and after calib -------------------"<<endl;
    cout << "Before calib: mean = " << h_before->GetMean() 
        << ", sigma = " << h_before->GetRMS() << endl;
    cout << "After calib:  mean = " << h_after->GetMean() 
        << ", sigma = " << h_after->GetRMS() << endl;
    // draw
    TCanvas *c = new TCanvas("c","Before vs After",800,600);
    h_before->SetLineColor(kRed);
    h_after ->SetLineColor(kBlue);
    h_before->Draw();
    h_after ->Draw("SAME");
    TLegend *leg=new TLegend(0.6,0.7,0.9,0.9);
    leg->AddEntry(h_before,"Before","l");
    leg->AddEntry(h_after, "After" ,"l");
    leg->Draw();
}
