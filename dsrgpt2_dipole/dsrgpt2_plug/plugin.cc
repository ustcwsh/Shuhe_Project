/*
 * @BEGIN LICENSE
 *
 * dsrgpt2_plug by Psi4 Developer, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Psi4; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */
#include "psi4/libdpd/dpd.h"
#include "psi4/libmints/vector.h" // <- needed to access SharedVector
#include "psi4/libmints/wavefunction.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"
#include "psi4/libmints/dipole.h"


#include <math.h>

// This allows us to be lazy in getting the spaces in DPD calls
#define ID(x) ints.DPD_ID(x)

double e=2.718281828;
double ss=5000;

namespace psi {
namespace dsrgpt2_plug {


extern "C" int read_options(std::string name, Options& options) {
    if (name == "dsrgpt2_plug" || options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);
    }

    return true;
}

void build_AOdipole_ints(SharedWavefunction wfn, SharedMatrix Dp) {

    std::shared_ptr<BasisSet> basisset = wfn->basisset();
    std::shared_ptr<IntegralFactory> ints_fac = std::make_shared<IntegralFactory>(basisset);
    int nbf = basisset->nbf();

    std::vector<SharedMatrix> AOdipole_ints_;
    // AOdipole_ints_.clear();
    for (const std::string& direction : {"X", "Y", "Z"}) {
        std::string name = "AO Dipole " + direction;
        AOdipole_ints_.push_back(SharedMatrix(new Matrix(name, nbf, nbf)));
    }
    std::shared_ptr<OneBodyAOInt> aodOBI(ints_fac->ao_dipole());
    aodOBI->compute(AOdipole_ints_);
    Dp->copy(AOdipole_ints_[2]);
}


extern "C" SharedWavefunction dsrgpt2_plug(SharedWavefunction ref_wfn, Options& options) {
    /*
     * This plugin shows a simple way of obtaining MO basis integrals, directly
     * from a DPD buffer.  It is also possible to generate integrals with labels
     * (IWL) formatted files, but that's not shown here.
     */
    int print = options.get_int("PRINT");

    // Grab the global (default) PSIO object, for file I/O
    std::shared_ptr<PSIO> psio(_default_psio_lib_);

    // Have the reference (SCF) wavefunction, ref_wfn
    if (!ref_wfn)
        throw PSIEXCEPTION("SCF has not been run yet!");

    // Quickly check that there are no open shell orbitals here...
    int nirrep = ref_wfn->nirrep();

    // For now, we'll just transform for closed shells and generate all integrals.
    // For more elaborate use of the LibTrans object, check out the plugin_mp2
    // example in the test suite.
    std::vector<std::shared_ptr<MOSpace>> spaces;
    spaces.push_back(MOSpace::all);
    IntegralTransform ints(ref_wfn, spaces, IntegralTransform::Restricted);
    ints.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
    // Use the IntegralTransform object's DPD instance, for convenience
    dpd_set_default(ints.get_dpd_id());

    /*
     * Now, loop over the DPD buffer, printing the integrals
     */
    dpdbuf4 K;
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    // To only process the permutationally unique integrals, change the
    // ID("[A,A]") to ID("[A>=A]+")
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[A,A]"), ID("[A,A]"), ID("[A>=A]+"),
                           ID("[A>=A]+"), 0, "MO Ints (AA|AA)");

    // 1. Read and store the two-electron integrals in chemist notation (pq|rs)
    // allocate a vector of size nmo^4
    size_t nmo = ref_wfn->nmo();
    size_t nmo4 = nmo * nmo * nmo * nmo;
    std::vector<double> mo_ints(nmo4, 0.0);

    // function to address a four-dimensional tensor of dimension dim * dim * dim * dim
    auto four_idx = [&](size_t p, size_t q, size_t r, size_t s, size_t dim) -> size_t {
        size_t dim2 = dim * dim;
        size_t dim3 = dim2 * dim;
        return (p * dim3 + q * dim2 + r * dim + s);
    };

    // read the integrals
    for (int h = 0; h < nirrep; ++h) {
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for (int pq = 0; pq < K.params->rowtot[h]; ++pq) {
            int p = K.params->roworb[h][pq][0];
            int q = K.params->roworb[h][pq][1];
            int psym = K.params->psym[p];
            int qsym = K.params->qsym[q];
            int prel = p - K.params->poff[psym];
            int qrel = q - K.params->qoff[qsym];
            for (int rs = 0; rs < K.params->coltot[h]; ++rs) {
                int r = K.params->colorb[h][rs][0];
                int s = K.params->colorb[h][rs][1];
                int rsym = K.params->rsym[r];
                int ssym = K.params->ssym[s];
                int rrel = r - K.params->roff[rsym];
                int srel = s - K.params->soff[ssym];
                // store the integrals
                mo_ints[four_idx(p, q, r, s, nmo)] = K.matrix[h][pq][rs];
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);
    psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    // 2. Build the antisymmetrized two-electron integrals in a spin orbital basis
    size_t nso = 2 * nmo;

    // define the order of spin orbitals and store it as a vector of pairs (orbital index,spin)
    std::vector<std::pair<size_t, int>> so_labels(nso);
    for (size_t n = 0; n < nmo; n++) {
        so_labels[2 * n] = std::make_pair(n, 0);     // 0 = alpha
        so_labels[2 * n + 1] = std::make_pair(n, 1); // 1 = beta
    }

    // allocate the vector that will store the spin orbital integrals
    size_t nso4 = nso * nso * nso * nso;
    std::vector<double> so_ints(nso4, 0.0);

    // form the integrals <pq||rs> = <pq|rs> - <pq|sr> = (pr|qs) - (ps|qr)
    for (size_t p = 0; p < nso; p++) {
        size_t p_orb = so_labels[p].first;
        int p_spin = so_labels[p].second;
        for (size_t q = 0; q < nso; q++) {
            size_t q_orb = so_labels[q].first;
            int q_spin = so_labels[q].second;
            for (size_t r = 0; r < nso; r++) {
                size_t r_orb = so_labels[r].first;
                int r_spin = so_labels[r].second;
                for (size_t s = 0; s < nso; s++) {
                    size_t s_orb = so_labels[s].first;
                    int s_spin = so_labels[s].second;

                    double integral = 0.0;
                    if ((p_spin == r_spin) and (q_spin == s_spin)) {
                        integral += mo_ints[four_idx(p_orb, r_orb, q_orb, s_orb, nmo)];
                    }
                    if ((p_spin == s_spin) and (q_spin == r_spin)) {
                        integral -= mo_ints[four_idx(p_orb, s_orb, q_orb, r_orb, nmo)];
                    }
                    so_ints[four_idx(p, q, r, s, nso)] = integral;
                }
            }
        }
    }

    // 3. Get the orbital energies from the reference wave function
    SharedVector epsilon_a = ref_wfn->epsilon_a();
    SharedVector epsilon_b = ref_wfn->epsilon_b();
    std::vector<double> epsilon(nso, 0.0);
    for (size_t p = 0; p < nso; p++) {
        size_t p_orb = so_labels[p].first;
        size_t p_spin = so_labels[p].second;
        if (p_spin == 0){
            epsilon[p] = epsilon_a->get(p_orb);
        }else{
            epsilon[p] = epsilon_b->get(p_orb);
        }
    }

    // 4. Form list of occupied and virtual orbitals
    int na = ref_wfn->nalpha();
    int nb = ref_wfn->nbeta();
    int nocc = na + nb;
    // ASSUMES RESTRICTED ORBITALS

    std::vector<size_t> O;
    std::vector<size_t> V;

    for (int i = 0; i < nocc; i++) {
        O.push_back(i);
    }
    for (int a = nocc; a < nso; a++) {
        V.push_back(a);
    }

    double mp2_energy = 0.0;
    for (int i : O) {
        for (int j : O) {
            for (int a : V) {
                for (int b : V) {
                    double Vijab = so_ints[four_idx(i, j, a, b, nso)];
                    double Dijab = epsilon[i] + epsilon[j] - epsilon[a] - epsilon[b];
                    mp2_energy += 0.25 * Vijab * Vijab / Dijab * ( 1 - pow(e,(-2.0 * ss * Dijab * Dijab))) ;
                }
            }
        }
    }

    std::cout << "what";

    double rhf_energy = ref_wfn->reference_energy();

    outfile->Printf("\n\n    ==> Spin orbital SR_DSRG_PT2 energy <==\n");
    outfile->Printf("    RHF total energy              %20.12f\n", rhf_energy);
    outfile->Printf("    SR_PT2 correlation energy     %20.12f\n", mp2_energy);
    outfile->Printf("    SR_DSRG_PT2 total energy      %20.12f\n", rhf_energy + mp2_energy);

    Process::environment.globals["CURRENT ENERGY"] = rhf_energy + mp2_energy;

//***************************************************
  
//     SharedMatrix Dp (new Matrix("Dipole correction matrix", 1, 7, 7, 0));  
//     build_AOdipole_ints(ref_wfn, Dp);
//     double dipole_energy = 0.0;

//     for(int i : O){
//         dipole_energy += Dp->get(0,i,i);
//     }

//     for (int i : O) {
//         for (int j : O) {
//             for (int a : V) {
//                 for (int b : V) {
//                     double Vijab = so_ints[four_idx(i, j, a, b, nso)];
//                     double Muijab = Dp->get(0,a,a) + Dp->get(0,b,b) - Dp->get(0,i,i) - Dp->get(0,j,j);
//                     dipole_energy += 0.25 * Vijab * Vijab / Muijab  ;
//                 }
//             }
//         }
//     }

// std::cout << dipole_energy;
std::cout << "what";





//***************************************************

    return ref_wfn;
}

} // namespace dsrgpt2_plug
} // namespace psi
