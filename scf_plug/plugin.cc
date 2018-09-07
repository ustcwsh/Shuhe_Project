#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/dipole.h"
#include "psi4/libpsio/psio.hpp"
#include <math.h>
#include <iomanip>
#include <vector>

double e = 2.718281828;

namespace psi{ namespace scf_plug {

extern "C" PSI_API
int read_options(std::string name, Options& options)
{
    if (name == "SCF_PLUG"|| options.read_globals()) {
        options.add_int("PRINT", 0);
        options.add_int("GRADIENT", 0);
        options.add_double("CVG", 0);
        options.add_double("PERT", 0);
        options.add_double("S", 0);
    }
    return true;
}

void FormDensityMatrix(SharedMatrix D, SharedMatrix C, int nmo, int occ)
{
    D->zero();

    for(int i = 0; i < nmo; ++i){
		for(int j = 0; j < nmo; ++j){
			for(int k = 0; k < occ; ++k){
				D->add(0, i, j, C->get(0, i, k) * C->get(0, j, k));
			}
		}
	}
}

void FormNewFockMatrix(SharedMatrix F, SharedMatrix H, SharedMatrix D, SharedMatrix eri, int nmo)
{
    for(int p = 0; p < nmo; ++p){
    	for(int q = 0; q < nmo; ++q){
    		F->set(0, p, q, H->get(0, p, q));
    		for(int r = 0; r < nmo; ++r){
    			for(int s = 0; s < nmo; ++s){
    				F->add(0, p, q, D->get(0, r, s) * ( 2.0 * eri->get(0, p * nmo + q, r * nmo + s) - eri->get(0, p * nmo + r, q * nmo + s)));
    			}
    		}
    	}
    }
}

double ElecEnergy(double Elec, SharedMatrix D, SharedMatrix H, SharedMatrix F, int nmo)
{
    Elec = 0;

    for(int p = 0; p < nmo; ++p){
    	for(int q = 0; q < nmo; ++q){
    		Elec += D->get(0, p, q)*(H->get(0, p, q)+F->get(0, p, q));
    	}
    }
    return Elec;
}

void AO2MO_TwoElecInts(SharedMatrix eri, SharedMatrix eri_mo, SharedMatrix C, int nmo)
{
    int dims[] = {0};
    dims[0] = nmo;
    SharedMatrix X (new Matrix("X", 1, dims, dims, 0));
    SharedMatrix Y (new Matrix("Y", 1, dims, dims, 0));
    SharedMatrix eri_temp = eri->clone();
                 eri_temp->zero();

    for(int i = 0; i < nmo; ++i)
    {
        for(int j = 0; j < nmo; ++j)
        {
            for(int k = 0; k < nmo; ++k)
            {
                for(int l = 0; l <= k; ++l)
                {
                    X->set(0, k, l, eri->get(0, i * nmo + j, k * nmo +l));
                    X->set(0, l, k, X->get(0, k, l));                 
                }
            }
            Y = Matrix::triplet(C, X, C, true, false, false);
            for(int k = 0; k < nmo; ++k)
            {
                for(int l = 0; l < nmo; ++l)
                {
                    eri_temp->set(0, k * nmo + l, i * nmo + j, Y->get(0, k, l));      
                }
            }        
        }
    }
    for(int k = 0; k < nmo; ++k)
    {
        for(int l = 0; l < nmo; ++l)
        {
            X->zero();
            Y->zero();
            for(int i = 0; i < nmo; ++i)
            {
                for(int j = 0; j <= i; ++j)
                {
                    X->set(0, i, j, eri_temp->get(0, k * nmo + l, i * nmo +j));
                    X->set(0, j, i, X->get(0, i, j));                     
                }
            }
            Y = Matrix::triplet(C, X, C, true, false, false);
            for(int i = 0; i < nmo; ++i)
            {
                for(int j = 0; j < nmo; ++j)
                {
                    eri_mo->set(0, k * nmo + l, i * nmo +j, Y->get(0, i, j));
                }
            }
        }
    }
}

void AO2MO_FockMatrix(SharedMatrix F, SharedMatrix F_MO, SharedMatrix C, int nmo)
{
    for(int i = 0; i < nmo; ++i)
    {
        for(int j = 0; j < nmo; ++j)
        {
            double ints_temp = 0.0;

            for(int p = 0; p < nmo; ++p)
            {
                for(int q = 0; q < nmo; ++q)
                {
                    ints_temp += F->get(0, p, q) * C->get(0, p, i) * C->get(0, q, j);
                }
            }
            F_MO->set(0, i, j, ints_temp);
        }
    }
}

double MP2_Energy_SO(SharedMatrix eri_mo, SharedMatrix F_MO, int nso, int doccpi, std::vector<double> so_ints, std::vector<double> epsilon_ijab, int frozen_c, int frozen_v)
{
    double Emp2 = 0.0;
    int idx;

    for(int i = frozen_c; i < 2 * doccpi; ++i)
    {
        for(int j = frozen_c; j < 2 * doccpi; ++j)
        {
            for(int a = 2 * doccpi; a < nso - frozen_v; ++a)
            {
                for(int b = 2 * doccpi; b < nso - frozen_v; ++b)
                {
                    idx = i * nso * nso * nso + j * nso * nso + a * nso + b;
                    Emp2 += 0.25 * so_ints[idx] * so_ints[idx] / epsilon_ijab[idx];
                }
            }
        }
    }
    return(Emp2);
}


/******************** TEST delete by Sep.1. ********************/

double MP2_Energy_MO(SharedMatrix eri_mo, SharedMatrix F_MO, int nmo, int doccpi, std::vector<double> mo_ints_aa, std::vector<double> mo_ints_bb, std::vector<double> mo_ints_ab, std::vector<double> epsilon_ijab_aa, int frozen_c, int frozen_v)
{
    double Emp2 = 0.0;
    int idx,idx1,idx2,idx3;

    for(int i = frozen_c/2; i <  doccpi; ++i)
    {
        for(int j = frozen_c/2; j <  doccpi; ++j)
        {
            for(int a = doccpi; a < nmo - frozen_v/2; ++a)
            {
                for(int b = doccpi; b < nmo - frozen_v/2; ++b)
                {
                    idx = i * nmo * nmo * nmo + j * nmo * nmo + a * nmo + b;
                    idx1 = i * nmo * nmo * nmo + j * nmo * nmo + b * nmo + a;
                    idx2 = j * nmo * nmo * nmo + i * nmo * nmo + a * nmo + b;
                    idx3 = j * nmo * nmo * nmo + i * nmo * nmo + b * nmo + a;
                    
                    Emp2 += 0.25 * mo_ints_aa[idx] * mo_ints_aa[idx] / epsilon_ijab_aa[idx];
                    Emp2 += 0.25 * mo_ints_bb[idx] * mo_ints_bb[idx] / epsilon_ijab_aa[idx];
                    Emp2 += 0.25 * mo_ints_ab[idx] * mo_ints_ab[idx] / epsilon_ijab_aa[idx];
                    Emp2 += 0.25 * mo_ints_ab[idx1] * mo_ints_ab[idx1] / epsilon_ijab_aa[idx];
                    Emp2 += 0.25 * mo_ints_ab[idx2] * mo_ints_ab[idx2] / epsilon_ijab_aa[idx];
                    Emp2 += 0.25 * mo_ints_ab[idx3] * mo_ints_ab[idx3] / epsilon_ijab_aa[idx];
                }
            }
        }
    }
    return(Emp2);
}
/******************** TEST delete by Sep.1. ********************/



double DSRG_PT2_Energy_SO(SharedMatrix eri_mo, SharedMatrix F_MO, int nso, int doccpi, std::vector<double> so_ints, std::vector<double> epsilon_ijab, double S, int frozen_c, int frozen_v)
{
    double Edsrgpt2 = 0.0;
    int idx;

    for(int i = frozen_c; i < 2 * doccpi; ++i)
    {
        for(int j = frozen_c; j < 2 * doccpi; ++j)
        {
            for(int a = 2 * doccpi; a < nso - frozen_v; ++a)
            {
                for(int b = 2 * doccpi; b < nso - frozen_v; ++b)
                {
                    idx = i * nso * nso * nso + j * nso * nso + a * nso + b;
                    Edsrgpt2 += 0.25 * so_ints[idx] * so_ints[idx] / epsilon_ijab[idx] * (1.0 - pow(e, -2.0 * S * epsilon_ijab[idx] * epsilon_ijab[idx]));
                }
            }
        }
    }
    return(Edsrgpt2);
}








/******************** TEST delete by Sep.2. ********************/



double DSRG_PT2_Energy_MO(SharedMatrix eri_mo, SharedMatrix F_MO, int nmo, int doccpi, std::vector<double> mo_ints_aa, std::vector<double> mo_ints_bb, std::vector<double> mo_ints_ab, std::vector<double> epsilon_ijab_aa, std::vector<double> epsilon_ijab_bb, std::vector<double> epsilon_ijab_ab, double S, int frozen_c, int frozen_v)
{
    double Edsrgpt2 = 0.0;
    int idx, idx1, idx2, idx3;

    for(int i = frozen_c/2; i < doccpi; ++i)
    {
        for(int j = frozen_c/2; j < doccpi; ++j)
        {
            for(int a = doccpi; a < nmo - frozen_v/2; ++a)
            {
                for(int b = doccpi; b < nmo - frozen_v/2; ++b)
                {

                    idx = i * nmo * nmo * nmo + j * nmo * nmo + a * nmo + b;
                    idx1 = i * nmo * nmo * nmo + j * nmo * nmo + b * nmo + a;
                    idx2 = j * nmo * nmo * nmo + i * nmo * nmo + a * nmo + b;
                    idx3 = j * nmo * nmo * nmo + i * nmo * nmo + b * nmo + a;
                    
                    Edsrgpt2 += 0.25 * mo_ints_aa[idx] * mo_ints_aa[idx] / epsilon_ijab_aa[idx] * (1.0 - pow(e, -2.0 * S * epsilon_ijab_aa[idx] * epsilon_ijab_aa[idx]));
                    Edsrgpt2 += 0.25 * mo_ints_bb[idx] * mo_ints_bb[idx] / epsilon_ijab_bb[idx] * (1.0 - pow(e, -2.0 * S * epsilon_ijab_bb[idx] * epsilon_ijab_bb[idx]));
                    Edsrgpt2 +=  1.0 * mo_ints_ab[idx] * mo_ints_ab[idx] / epsilon_ijab_ab[idx] * (1.0 - pow(e, -2.0 * S * epsilon_ijab_ab[idx] * epsilon_ijab_ab[idx]));
                    // Edsrgpt2 += 0.25 * mo_ints_ab[idx1] * mo_ints_ab[idx1] / epsilon_ijab_ab[idx] * (1.0 - pow(e, -2.0 * S * epsilon_ijab_ab[idx] * epsilon_ijab_ab[idx]));
                    // Edsrgpt2 += 0.25 * mo_ints_ab[idx2] * mo_ints_ab[idx2] / epsilon_ijab_ab[idx] * (1.0 - pow(e, -2.0 * S * epsilon_ijab_ab[idx] * epsilon_ijab_ab[idx]));
                    // Edsrgpt2 += 0.25 * mo_ints_ab[idx3] * mo_ints_ab[idx3] / epsilon_ijab_ab[idx] * (1.0 - pow(e, -2.0 * S * epsilon_ijab_ab[idx] * epsilon_ijab_ab[idx]));
                


                }
            }
        }
    }
    return(Edsrgpt2);
}
/******************** TEST delete by Sep.2. ********************/


















void build_AOdipole_ints(SharedWavefunction wfn, SharedMatrix Dp) 
{
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

extern "C" PSI_API
SharedWavefunction scf_plug(SharedWavefunction ref_wfn, Options& options)
{
    //Parameters Declaration
    std::shared_ptr<Molecule> molecule = Process::environment.molecule();
    molecule->update_geometry();
    molecule->print();
    std::shared_ptr<BasisSet> ao_basisset = ref_wfn->basisset();
    MintsHelper mints(ref_wfn);
    int      dims[] = {ao_basisset->nbf()};
    int      print = options.get_int("PRINT");
    int      gradient = options.get_int("GRADIENT");
    double   CVG = options.get_double("CVG");
    double   pert = options.get_double("PERT");
    double   S_const = options.get_double("S");
    int      iternum = 1;
    double   energy_pre;
    int      doccpi = 0;
    double   Enuc = molecule->nuclear_repulsion_energy(ref_wfn->get_dipole_field_strength());
    double   Elec, Etot, Emp2, Edsrg_pt2;
    int      irrep_num = ref_wfn->nirrep();
    int      nmo = dims[0];
    size_t   nso = 2 * dims[0];
    int      frozen_c = 2;
    int      frozen_v = 14;
    std::shared_ptr<MatrixFactory> factory(new MatrixFactory);
    factory->init_with(1, dims, dims);
    SharedMatrix overlap = mints.ao_overlap();
    SharedMatrix kinetic = mints.ao_kinetic();
    SharedMatrix potential = mints.ao_potential();

    SharedMatrix F_a = ref_wfn->Fa();
    SharedMatrix F_b = ref_wfn->Fb();
    SharedMatrix C_a = ref_wfn->Ca();
    SharedMatrix C_b = ref_wfn->Cb();
    SharedMatrix F_MO_a (new Matrix("Fock_MO_alpha matrix", 1, dims, dims, 0));
    SharedMatrix F_MO_b (new Matrix("Fock_MO_beta matrix", 1, dims, dims, 0));




	SharedMatrix Omega = factory->create_shared_matrix("Omega");
    SharedMatrix F (new Matrix("Fock matrix", 1, dims, dims, 0));
    SharedMatrix F_uptp (new Matrix("Unperturbed Fock matrix", 1, dims, dims, 0));
    SharedMatrix F_MO (new Matrix("Fock_MO matrix", 1, dims, dims, 0));
    SharedMatrix F_MO_uptp (new Matrix("unperturbed Fock_MO matrix", 1, dims, dims, 0));
    SharedMatrix C (new Matrix("C matrix", 1, dims, dims, 0));
    SharedMatrix C_uptp (new Matrix("Unperturbed C matrix", 1, dims, dims, 0));
    SharedMatrix D (new Matrix("Density matrix", 1, dims, dims, 0));
    SharedMatrix D_uptp (new Matrix("Unperturbed Density matrix", 1, dims, dims, 0));
    SharedMatrix S (new Matrix("S matrix", 1, dims, dims, 0));
    SharedMatrix H = factory->create_shared_matrix("H");
    SharedMatrix H_uptb = factory->create_shared_matrix("Unperturbed H");
    SharedMatrix eri = mints.ao_eri();
    SharedMatrix eri_mo = eri->clone();
                 eri_mo->zero();
    SharedMatrix Dp (new Matrix("Dipole correction matrix", 1, dims, dims, 0));
    SharedMatrix Dp_temp (new Matrix("Dipole correction matrix Copy", 1, dims, dims, 0));
    SharedMatrix evecs (new Matrix("evecs", 1, dims, dims, 0));
    SharedVector evals (new Vector("evals", 1, dims));
    SharedVector ndip = DipoleInt::nuclear_contribution(Process::environment.molecule(), Vector3(0.0, 0.0, 0.0));

        
    Enuc += pert * ndip->get(0, 2);
    Dimension doccpi_add = ref_wfn->doccpi();

    for(int i = 0; i < irrep_num; ++i)
    {
    	doccpi += doccpi_add[i];
    }

    std::cout << std::endl;

/************************ index ************************/

    auto four_idx = [&](size_t p, size_t q, size_t r, size_t s, size_t dim) -> size_t 
    {
        size_t dim2 = dim * dim;
        size_t dim3 = dim2 * dim;
        return (p * dim3 + q * dim2 + r * dim + s);
    };

    auto two_idx = [&](size_t p, size_t q, size_t dim) -> size_t 
    {
        return (p * dim + q);
    };

/************************ SCF ************************/

    //Create H matrix
    H->copy(kinetic);
    H->add(potential);
    H_uptb->copy(H);
    build_AOdipole_ints(ref_wfn, Dp);
    Dp_temp->copy(Dp);
    Dp->Matrix::scale(pert);
    H->add(Dp);
    Dp->copy(Dp_temp);

    //Create S^(-1/2) Matrix
    Omega->zero();
    overlap->diagonalize(evecs, evals);

    for(int i = 0; i < nmo; ++i)
    {
    	Omega->set(0, i, i, 1.0 / sqrt(evals->get(0, i)));
    }
     
    S = Matrix::triplet(evecs, Omega, evecs, false, false, true);

    //Create original fock matrix using transformation on H
	F = Matrix::triplet(S, H, S, true, false, false);
    F_uptp = Matrix::triplet(S, H_uptb, S, true, false, false);

    //Create C matrix
    F->diagonalize(evecs, evals);
    C = Matrix::doublet(S, evecs, false, false);
    F_uptp->diagonalize(evecs, evals);
    C_uptp = Matrix::doublet(S, evecs, false, false);    

    //Create Density matrix
    FormDensityMatrix(D, C, nmo, doccpi);
    FormDensityMatrix(D_uptp, C_uptp, nmo, doccpi);

    //Create new Fock matrix
	FormNewFockMatrix(F, H, D, eri, nmo);
    FormNewFockMatrix(F_uptp, H_uptb, D_uptp, eri, nmo);

    //Calculate the energy
	Elec = ElecEnergy(Elec, D_uptp, H_uptb, F_uptp, nmo);/*!!!!! TEST !!!! unperturbed*/
    Etot = Elec + Enuc;
    energy_pre = Etot; 
    Etot += CVG + 1.0;

    //SCF iteration

// while( fabs(energy_pre - Etot) > CVG){
    iternum = 0;

    while(iternum < 100)
    {
	    energy_pre = Etot;
        F = Matrix::triplet(S, F, S, true, false, false);
        F->diagonalize(evecs, evals);
        C = Matrix::doublet(S, evecs, false, false);
        FormDensityMatrix(D, C, nmo, doccpi);
	    FormNewFockMatrix(F, H, D, eri, nmo);
        iternum++;
        /*********** unperturbed C *********/
        F_uptp = Matrix::triplet(S, F_uptp, S, true, false, false);
        F_uptp->diagonalize(evecs, evals);
        C_uptp = Matrix::doublet(S, evecs, false, false);
        FormDensityMatrix(D_uptp, C_uptp, nmo, doccpi);
        FormNewFockMatrix(F_uptp, H_uptb, D_uptp, eri, nmo);
        Elec = ElecEnergy(Elec, D_uptp, H_uptb, F_uptp, nmo);
        Etot = Elec + Enuc;
    }

    double Escf = Etot;


/************************ MP2 & DSRG-PT2 (Orbital irrelevant) ver1.0 ************************/
   
    SharedMatrix Dp_d (new Matrix("dipole diagonal matrix", 1, dims, dims, 0));
    Dp_d->zero();
    SharedMatrix Dp_mo = Dp->clone();    
    AO2MO_FockMatrix(Dp, Dp_mo, C_uptp, nmo);    
    /* dipole OV,OO,VV block */
    Dp_d->copy(Dp_mo);
    Dp_d->Matrix::scale(pert);
    AO2MO_FockMatrix(F_uptp, F_MO, C_uptp, nmo);
    F_MO_uptp->copy(F_MO);
    F_MO->add(Dp_d);



     /************** TEST ****************/

    AO2MO_FockMatrix(F_a, F_MO_a, C_a, nmo);
    AO2MO_FockMatrix(F_b, F_MO_b, C_b, nmo);

     /************** TEST ****************/
    // F_MO_b->print();
    // F_MO_a->print();





/************************ MP2 & DSRG-PT2 (Orbital irrelevant) ver3.0 O-O block ************************/

    //dipole OO,VV,OV block 
    F_MO->diagonalize(evecs, evals);
    SharedMatrix F_MO_1 (new Matrix("rotated Fock matrix", 1, dims, dims, 0));
    F_MO_1->zero();

    for(int i = 0; i < nmo; ++i){
        F_MO_1->set(0, i, i, evals->get(0, i));
    }

    SharedMatrix C_pert (new Matrix("C_pert matrix", 1, dims, dims, 0));
    SharedMatrix U_pert (new Matrix("transformation matrix", 1, dims, dims, 0));

    U_pert->copy(evecs);
    C_pert = Matrix::doublet(C_uptp, U_pert, false, false);

/********** Obtain the gradient or not  *************/

    if(gradient)
    {
        F_MO->copy(F_MO_uptp);
    }
    else
    {
        C_uptp->copy(C);
        AO2MO_FockMatrix(F, F_MO, C, nmo);    
    }

/************************ MP2 & DSRG-PT2 (Orbital irrelevant) ver2.0 ************************/

    AO2MO_TwoElecInts(eri, eri_mo, C_uptp, nmo);

    // define the order of spin orbitals and store it as a vector of pairs (orbital index,spin)
    std::vector<std::pair<size_t, int>> so_labels(nso);
    for (size_t n = 0; n < nmo; n++) {
        so_labels[2 * n] = std::make_pair(n, 0);     // 0 = alpha
        so_labels[2 * n + 1] = std::make_pair(n, 1); // 1 = beta
    }

    // allocate the vector that will store the spin orbital integrals
    size_t nso2 = nso * nso;
    size_t nso4 = nso2 * nso2;
    std::vector<double> so_ints(nso4, 0.0);
    std::vector<double> amp_t(nso4, 0.0);
    std::vector<double> amp_t_dsrg(nso4, 0.0);
    std::vector<double> epsilon(nso, 0.0);
    std::vector<double> epsilon_ijab(nso4, 0.0);


    for (size_t p = 0; p < nso; ++p){
        epsilon[p] = F_MO->get(0, p / 2, p / 2);
    }

    for (size_t i = 0; i < 2 * doccpi; ++i)
    {
        for (size_t j = 0; j < 2 * doccpi; ++j)
        {
            for (size_t a = 2 * doccpi; a < nso; ++a)
            {
                for (size_t b = 2 * doccpi; b < nso; ++b)
                {
                    epsilon_ijab[four_idx(i, j, a, b, nso)] = epsilon[i] + epsilon[j] - epsilon[a] - epsilon[b];
                }
            }
        }
    }



    /************  TEST  delete by Sept.1. ************/

    size_t nmo2 = nmo * nmo;
    size_t nmo4 = nmo2 * nmo2;
    std::vector<double> mo_ints_aa(nmo4, 0.0);
    std::vector<double> mo_ints_bb(nmo4, 0.0);
    std::vector<double> mo_ints_ab(nmo4, 0.0);   // V_{abab}

    std::vector<double> epsilon_a(nmo, 0.0);
    std::vector<double> epsilon_ijab_aa(nmo4, 0.0);

/*********** Sept 6 ****************/
    std::vector<double> epsilon_b(nmo, 0.0);
    std::vector<double> epsilon_ijab_bb(nmo4, 0.0);
     std::vector<double> epsilon_ijab_ab(nmo4, 0.0);


    for (size_t p = 0; p < nmo; ++p){
        epsilon_a[p] = F_MO_a->get(0, p, p);
        epsilon_b[p] = F_MO_b->get(0, p, p);
    }
    

    for (size_t i = 0; i < doccpi; ++i)
    {
        for (size_t j = 0; j < doccpi; ++j)
        {
            for (size_t a = doccpi; a < nmo; ++a)
            {
                for (size_t b = doccpi; b < nmo; ++b)
                {
                    epsilon_ijab_aa[four_idx(i, j, a, b, nmo)] = epsilon_a[i] + epsilon_a[j] - epsilon_a[a] - epsilon_a[b];
                    epsilon_ijab_bb[four_idx(i, j, a, b, nmo)] = epsilon_b[i] + epsilon_b[j] - epsilon_b[a] - epsilon_b[b];
                    epsilon_ijab_ab[four_idx(i, j, a, b, nmo)] = epsilon_a[i] + epsilon_b[j] - epsilon_a[a] - epsilon_b[b];
                }
            }
        }
    }


    // for (size_t i = 0; i < doccpi; ++i)
    // {
    //     for (size_t j = 0; j < doccpi; ++j)
    //     {
    //         for (size_t a = doccpi; a < nmo; ++a)
    //         {
    //             for (size_t b = doccpi; b < nmo; ++b)
    //             {
    //                 epsilon_bb[four_idx(i, j, a, b, nmo)] = epsilon[2 * i + 1] + epsilon[2 * j + 1] - epsilon[2 * a + 1] - epsilon[2 * b + 1];
    //             }
    //         }
    //     }
    // }

    // for (size_t i = 0; i < doccpi; ++i)
    // {
    //     for (size_t j = 0; j < doccpi; ++j)
    //     {
    //         for (size_t a = doccpi; a < nmo; ++a)
    //         {
    //             for (size_t b = doccpi; b < nmo; ++b)
    //             {
    //                 epsilon_ab[four_idx(i, j, a, b, nmo)] = epsilon[2 * i] + epsilon[2 * j + 1] - epsilon[2 * a] - epsilon[2 * b + 1];
    //             }
    //         }
    //     }
    // }



    std::vector<double> amp_t_dsrg_aa(nmo4, 0.0);
    std::vector<double> amp_t_dsrg_bb(nmo4, 0.0);
    std::vector<double> amp_t_dsrg_ab(nmo4, 0.0);


    //form the integrals <pq||rs> = <pq|rs> - <pq|sr> = (pr|qs) - (ps|qr)
    for (size_t p = 0; p < nmo; p++) 
    {
        for (size_t q = 0; q < nmo; q++) 
        {
            for (size_t r = 0; r < nmo; r++) 
            {
                for (size_t s = 0; s < nmo; s++) 
                {
                        mo_ints_aa[four_idx(p, q, r, s, nmo)] = eri_mo->get(0, p * nmo + r, q * nmo + s) - eri_mo->get(0, p * nmo + s, q * nmo + r);

                        mo_ints_bb[four_idx(p, q, r, s, nmo)] = eri_mo->get(0, p * nmo + r, q * nmo + s) - eri_mo->get(0, p * nmo + s, q * nmo + r);
        
                        mo_ints_ab[four_idx(p, q, r, s, nmo)] = eri_mo->get(0, p * nmo + r, q * nmo + s);
 
                    if(p < doccpi && q < doccpi && r >= doccpi && s >= doccpi)
                    {
                        amp_t_dsrg_aa[four_idx(p, q, r, s, nmo)] = mo_ints_aa[four_idx(p, q, r, s, nmo)] / epsilon_ijab_aa[four_idx(p, q, r, s, nmo)] * (1.0 - pow(e, -S_const * epsilon_ijab_aa[four_idx(p, q, r, s, nmo)] * epsilon_ijab_aa[four_idx(p, q, r, s, nmo)]));
                        amp_t_dsrg_bb[four_idx(p, q, r, s, nmo)] = mo_ints_bb[four_idx(p, q, r, s, nmo)] / epsilon_ijab_bb[four_idx(p, q, r, s, nmo)] * (1.0 - pow(e, -S_const * epsilon_ijab_bb[four_idx(p, q, r, s, nmo)] * epsilon_ijab_bb[four_idx(p, q, r, s, nmo)]));
                        amp_t_dsrg_ab[four_idx(p, q, r, s, nmo)] = mo_ints_ab[four_idx(p, q, r, s, nmo)] / epsilon_ijab_ab[four_idx(p, q, r, s, nmo)] * (1.0 - pow(e, -S_const * epsilon_ijab_ab[four_idx(p, q, r, s, nmo)] * epsilon_ijab_ab[four_idx(p, q, r, s, nmo)]));
                    }
                }
            }
        }
    }










    /************  TEST  delete by Sept.1. ************/







    //form the integrals <pq||rs> = <pq|rs> - <pq|sr> = (pr|qs) - (ps|qr)
    for (size_t p = 0; p < nso; p++) 
    {
        size_t p_orb = so_labels[p].first;
        int p_spin = so_labels[p].second;
        for (size_t q = 0; q < nso; q++) 
        {
            size_t q_orb = so_labels[q].first;
            int q_spin = so_labels[q].second;
            for (size_t r = 0; r < nso; r++) 
            {
                size_t r_orb = so_labels[r].first;
                int r_spin = so_labels[r].second;
                for (size_t s = 0; s < nso; s++) 
                {
                    size_t s_orb = so_labels[s].first;
                    int s_spin = so_labels[s].second;
                    double integral = 0.0;

                    if ((p_spin == r_spin) and (q_spin == s_spin)) 
                    {
                        integral += eri_mo->get(0, p / 2 * nmo + r / 2, q / 2 * nmo + s / 2);
                    }    
                    if ((p_spin == s_spin) and (q_spin == r_spin)) 
                    {
                        integral -= eri_mo->get(0, p / 2 * nmo + s / 2, q / 2 * nmo + r / 2);
                    }
                    so_ints[four_idx(p, q, r, s, nso)] = integral;
                    if(p < 2 * doccpi && q < 2 * doccpi && r >= 2 * doccpi && s >= 2 * doccpi)
                    {
                        amp_t[four_idx(p, q, r, s, nso)] = integral / epsilon_ijab[four_idx(p, q, r, s, nso)];
                        amp_t_dsrg[four_idx(p, q, r, s, nso)] = integral / epsilon_ijab[four_idx(p, q, r, s, nso)] * (1.0 - pow(e, -S_const * epsilon_ijab[four_idx(p, q, r, s, nso)] * epsilon_ijab[four_idx(p, q, r, s, nso)]));
                    }
                }
            }
        }
    }


    int dims_nso2[] = {0};
    dims_nso2[0] = nso;
    // SharedMatrix D_MP2 (new Matrix("MP2 Dipole Density matrix", 1, dims_nso2, dims_nso2, 0));
    SharedMatrix Z_MP2 (new Matrix("Z MP2 matrix", 1, dims_nso2, dims_nso2, 0));
    SharedMatrix Z_temp (new Matrix("Z matrix temporal", 1, dims_nso2, dims_nso2, 0));
    SharedMatrix D_MP2 (new Matrix("MP2 Dipole Density matrix", 1, dims_nso2, dims_nso2, 0));

    Z_MP2->zero();
    D_MP2->zero();





    // for(int i = frozen_c; i < 2 * doccpi; ++i)
    // {
    //     for(int j = frozen_c; j < 2 * doccpi; ++j)
    //     {
    //         for(int a = 2 * doccpi; a < nso - frozen_v; ++a)
    //         {
    //             for(int b = 2 * doccpi; b < nso - frozen_v; ++b)
    //             {
    //                 double temp;

    //                 temp = - 0.5 * amp_t_dsrg[four_idx(i, j, a, b, nso)] * amp_t_dsrg[four_idx(i, j, a, b, nso)] * (1.0 + pow(e, -S_const * epsilon_ijab[four_idx(i, j, a, b, nso)] * epsilon_ijab[four_idx(i, j, a, b, nso)])) / (1.0 - pow(e, -S_const * epsilon_ijab[four_idx(i, j, a, b, nso)] * epsilon_ijab[four_idx(i, j, a, b, nso)])) + 2.0 * S_const * so_ints[four_idx(i, j, a, b, nso)] * so_ints[four_idx(i, j, a, b, nso)] * pow(e, -2.0 * S_const * epsilon_ijab[four_idx(i, j, a, b, nso)] * epsilon_ijab[four_idx(i, j, a, b, nso)]);
    //                 D_MP2->add(0, i, i, temp);
    //                 D_MP2->add(0, a, a, -temp);
    //             }
    //         }
    //     }
    // }
    /********************************** TEST (DONE) **********************************/
    for(int i = frozen_c/2; i < doccpi; ++i)
    {
        for(int j = frozen_c/2; j < doccpi; ++j)
        {
            for(int a = doccpi; a < nmo - frozen_v/2; ++a)
            {
                for(int b = doccpi; b < nmo - frozen_v/2; ++b)
                {
                    double temp1 ;
                    double temp2 ;
                    double temp3 ;

                    temp1 = -0.5 *  amp_t_dsrg_aa[four_idx(i, j, a, b, nmo)] * amp_t_dsrg_aa[four_idx(i, j, a, b, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_aa[four_idx(i, j, a, b, nmo)] * epsilon_ijab_aa[four_idx(i, j, a, b, nmo)])) / (1.0 - pow(e, -S_const * epsilon_ijab_aa[four_idx(i, j, a, b, nmo)] * epsilon_ijab_aa[four_idx(i, j, a, b, nmo)])) + 2.0 * S_const * mo_ints_aa[four_idx(i, j, a, b, nmo)] * mo_ints_aa[four_idx(i, j, a, b, nmo)] * pow(e, -2.0 * S_const * epsilon_ijab_aa[four_idx(i, j, a, b, nmo)] * epsilon_ijab_aa[four_idx(i, j, a, b, nmo)]);
                    temp2 = -0.5 *  amp_t_dsrg_bb[four_idx(i, j, a, b, nmo)] * amp_t_dsrg_bb[four_idx(i, j, a, b, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_bb[four_idx(i, j, a, b, nmo)] * epsilon_ijab_bb[four_idx(i, j, a, b, nmo)])) / (1.0 - pow(e, -S_const * epsilon_ijab_bb[four_idx(i, j, a, b, nmo)] * epsilon_ijab_bb[four_idx(i, j, a, b, nmo)])) + 2.0 * S_const * mo_ints_bb[four_idx(i, j, a, b, nmo)] * mo_ints_bb[four_idx(i, j, a, b, nmo)] * pow(e, -2.0 * S_const * epsilon_ijab_bb[four_idx(i, j, a, b, nmo)] * epsilon_ijab_bb[four_idx(i, j, a, b, nmo)]);
                    temp3 = 2.0 * ( -0.5 * amp_t_dsrg_ab[four_idx(i, j, a, b, nmo)] * amp_t_dsrg_ab[four_idx(i, j, a, b, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)] * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)])) / (1.0 - pow(e, -S_const * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)] * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)])) + 2.0 * S_const * mo_ints_ab[four_idx(i, j, a, b, nmo)] * mo_ints_ab[four_idx(i, j, a, b, nmo)] * pow(e, -2.0 * S_const * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)] * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)]));
                    D_MP2->add(0, 2*i, 2*i, temp1 + temp3);
                    D_MP2->add(0, 2*i+1, 2*i+1, temp2 + temp3);
                    D_MP2->add(0, 2*a, 2*a, -temp1 - temp3);
                    D_MP2->add(0, 2*a+1, 2*a+1, -temp2 - temp3);
                }
            }
        }
    }





    // for(int m = frozen_c; m < 2 * doccpi; ++m)
    // {
    //     for(int n = frozen_c; n < 2 * doccpi; ++n)
    //     {
    //         for(int j = frozen_c; j < 2 * doccpi; ++j)
    //         {
    //             for(int a = 2 * doccpi; a < nso - frozen_v; ++a)
    //             {
    //                 for(int b = 2 * doccpi; b < nso - frozen_v; ++b)
    //                 {
    //                     if(m!=n)
    //                     {
    //                         if(fabs(epsilon[m]-epsilon[n]) > 1e-6)
    //                         {
    //                             double t3;

    //                             t3 = 1.0 / (epsilon[n] - epsilon[m]) * so_ints[four_idx(m, j, a, b, nso)] * so_ints[four_idx(a, b, n, j, nso)] * ((1.0 - pow(e, -2.0 * S_const * epsilon_ijab[four_idx(n, j, a, b, nso)] * epsilon_ijab[four_idx(n, j, a, b, nso)])) / epsilon_ijab[four_idx(n, j, a, b, nso)] - (1.0 - pow(e, -2.0 * S_const * epsilon_ijab[four_idx(m, j, a, b, nso)] * epsilon_ijab[four_idx(m, j, a, b, nso)])) / epsilon_ijab[four_idx(m, j, a, b, nso)]);
    //                             Z_MP2->add(0, n, m, t3);
    //                         }
    //                         else
    //                         {
    //                             double t3;

    //                             t3 = so_ints[four_idx(m, j, a, b, nso)] * so_ints[four_idx(a, b, n, j, nso)] * (4.0 * S_const * pow(e, -2.0 * S_const * epsilon_ijab[four_idx(m, j, a, b, nso)] * epsilon_ijab[four_idx(m, j, a, b, nso)]) - (1.0 - pow(e, -2.0 * S_const * epsilon_ijab[four_idx(m, j, a, b, nso)] * epsilon_ijab[four_idx(m, j, a, b, nso)])) / epsilon_ijab[four_idx(m, j, a, b, nso)] / epsilon_ijab[four_idx(m, j, a, b, nso)]);
    //                             Z_MP2->add(0, n, m, t3);
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }     
    /********************************** TEST (DONE) **********************************/
    for(int m = frozen_c/2; m < doccpi; ++m)
    {
        for(int n = frozen_c/2; n < doccpi; ++n)
        {
            for(int j = frozen_c/2; j < doccpi; ++j)
            {
                for(int a = doccpi; a < nmo - frozen_v/2; ++a)
                {
                    for(int b = doccpi; b < nmo - frozen_v/2; ++b)
                    {
                        if(m != n)
                        {
                            if(fabs(epsilon_a[m]-epsilon_a[n]) > 1e-8)
                            {
                                double temp1;
                                double temp2;
                                double temp3;

                                temp1  = 1.0 / (epsilon_a[n] - epsilon_a[m]) * mo_ints_aa[four_idx(m, j, a, b, nmo)] * mo_ints_aa[four_idx(a, b, n, j, nmo)] * ((1.0 - pow(e, -2.0 * S_const * epsilon_ijab_aa[four_idx(n, j, a, b, nmo)] * epsilon_ijab_aa[four_idx(n, j, a, b, nmo)])) / epsilon_ijab_aa[four_idx(n, j, a, b, nmo)] - (1.0 - pow(e, -2.0 * S_const * epsilon_ijab_aa[four_idx(m, j, a, b, nmo)] * epsilon_ijab_aa[four_idx(m, j, a, b, nmo)])) / epsilon_ijab_aa[four_idx(m, j, a, b, nmo)]);
                                temp2  = 1.0 / (epsilon_a[n] - epsilon_a[m]) * mo_ints_bb[four_idx(m, j, a, b, nmo)] * mo_ints_bb[four_idx(a, b, n, j, nmo)] * ((1.0 - pow(e, -2.0 * S_const * epsilon_ijab_bb[four_idx(n, j, a, b, nmo)] * epsilon_ijab_bb[four_idx(n, j, a, b, nmo)])) / epsilon_ijab_bb[four_idx(n, j, a, b, nmo)] - (1.0 - pow(e, -2.0 * S_const * epsilon_ijab_bb[four_idx(m, j, a, b, nmo)] * epsilon_ijab_bb[four_idx(m, j, a, b, nmo)])) / epsilon_ijab_bb[four_idx(m, j, a, b, nmo)]);
                                temp3  = 2.0 / (epsilon_a[n] - epsilon_a[m]) * mo_ints_ab[four_idx(m, j, a, b, nmo)] * mo_ints_ab[four_idx(a, b, n, j, nmo)] * ((1.0 - pow(e, -2.0 * S_const * epsilon_ijab_ab[four_idx(n, j, a, b, nmo)] * epsilon_ijab_ab[four_idx(n, j, a, b, nmo)])) / epsilon_ijab_ab[four_idx(n, j, a, b, nmo)] - (1.0 - pow(e, -2.0 * S_const * epsilon_ijab_ab[four_idx(m, j, a, b, nmo)] * epsilon_ijab_ab[four_idx(m, j, a, b, nmo)])) / epsilon_ijab_ab[four_idx(m, j, a, b, nmo)]);
                                Z_MP2->add(0, 2*n, 2*m, temp1 + temp3);
                                Z_MP2->add(0, 2*n+1, 2*m+1, temp2 + temp3);
                            }
                            else
                            {
                                double temp1;
                                double temp2;
                                double temp3;

                                temp1 = mo_ints_aa[four_idx(m, j, a, b, nmo)] * mo_ints_aa[four_idx(a, b, n, j, nmo)] * (4.0 * S_const * pow(e, -2.0 * S_const * epsilon_ijab_aa[four_idx(m, j, a, b, nmo)] * epsilon_ijab_aa[four_idx(m, j, a, b, nmo)]) - (1.0 - pow(e, -2.0 * S_const * epsilon_ijab_aa[four_idx(m, j, a, b, nmo)] * epsilon_ijab_aa[four_idx(m, j, a, b, nmo)])) / epsilon_ijab_aa[four_idx(m, j, a, b, nmo)] / epsilon_ijab_aa[four_idx(m, j, a, b, nmo)]);
                                temp2 = mo_ints_bb[four_idx(m, j, a, b, nmo)] * mo_ints_bb[four_idx(a, b, n, j, nmo)] * (4.0 * S_const * pow(e, -2.0 * S_const * epsilon_ijab_bb[four_idx(m, j, a, b, nmo)] * epsilon_ijab_bb[four_idx(m, j, a, b, nmo)]) - (1.0 - pow(e, -2.0 * S_const * epsilon_ijab_bb[four_idx(m, j, a, b, nmo)] * epsilon_ijab_bb[four_idx(m, j, a, b, nmo)])) / epsilon_ijab_bb[four_idx(m, j, a, b, nmo)] / epsilon_ijab_bb[four_idx(m, j, a, b, nmo)]);
                                temp3 = mo_ints_ab[four_idx(m, j, a, b, nmo)] * mo_ints_ab[four_idx(a, b, n, j, nmo)] * (4.0 * S_const * pow(e, -2.0 * S_const * epsilon_ijab_ab[four_idx(m, j, a, b, nmo)] * epsilon_ijab_ab[four_idx(m, j, a, b, nmo)]) - (1.0 - pow(e, -2.0 * S_const * epsilon_ijab_ab[four_idx(m, j, a, b, nmo)] * epsilon_ijab_ab[four_idx(m, j, a, b, nmo)])) / epsilon_ijab_ab[four_idx(m, j, a, b, nmo)] / epsilon_ijab_ab[four_idx(m, j, a, b, nmo)]);
                                Z_MP2->add(0, 2*n, 2*m, temp1 + temp3);
                                Z_MP2->add(0, 2*n+1, 2*m+1, temp2 + temp3);
                            }
                        }
                    }
                }
            }
        }
    }






    // for(int c = 2 * doccpi; c < nso - frozen_v; ++c)
    // {
    //     for(int d = 2 * doccpi; d < nso - frozen_v; ++d)
    //     {
    //         for(int a = 2 * doccpi; a < nso - frozen_v; ++a)
    //         {
    //             for(int i = frozen_c; i < 2 * doccpi; ++i)
    //             {
    //                 for(int j = frozen_c; j < 2 * doccpi; ++j)
    //                 {
    //                     if(c != d)
    //                     {
    //                         if(fabs(epsilon[c] - epsilon[d]) > 1e-6)
    //                         {
    //                             double t4;

    //                             t4 = 1.0 / (epsilon[d] - epsilon[c]) * amp_t_dsrg[four_idx(i, j, a, c, nso)] * amp_t_dsrg[four_idx(i, j, a, d, nso)] * (epsilon_ijab[four_idx(i, j, a, c, nso)] * (1.0 + pow(e, -S_const * epsilon_ijab[four_idx(i, j, a, d, nso)] * epsilon_ijab[four_idx(i, j, a, d, nso)])) / (1.0 - pow(e, -S_const * epsilon_ijab[four_idx(i, j, a, c, nso)] * epsilon_ijab[four_idx(i, j, a, c, nso)])) - epsilon_ijab[four_idx(i, j, a, d, nso)] * (1.0 + pow(e, -S_const * epsilon_ijab[four_idx(i, j, a, c, nso)] * epsilon_ijab[four_idx(i, j, a, c, nso)])) / (1.0 - pow(e, -S_const * epsilon_ijab[four_idx(i, j, a, d, nso)] * epsilon_ijab[four_idx(i, j, a, d, nso)])));
    //                             Z_MP2->add(0, d, c, t4);
    //                         }
    //                         else
    //                         {
    //                             double t4;

    //                             t4 = so_ints[four_idx(i, j, a, c, nso)] * so_ints[four_idx(i, j, a, d, nso)] * (-4.0 * S_const * pow(e, -2.0 * S_const * epsilon_ijab[four_idx(i, j, a, c, nso)] * epsilon_ijab[four_idx(i, j, a, c, nso)]) + (1.0 - pow(e, -2.0 * S_const * epsilon_ijab[four_idx(i, j, a, c, nso)] * epsilon_ijab[four_idx(i, j, a, c, nso)])) / epsilon_ijab[four_idx(i, j, a, c, nso)] / epsilon_ijab[four_idx(i, j, a, c, nso)]);
    //                             Z_MP2->add(0, d, c, t4);
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }        
    /********************************** TEST (DONE) **********************************/
    for(int c = doccpi; c < nmo - frozen_v/2; ++c)
    {
        for(int d = doccpi; d < nmo - frozen_v/2; ++d)
        {
            for(int a = doccpi; a < nmo - frozen_v/2; ++a)
            {
                for(int i = frozen_c/2; i < doccpi; ++i)
                {
                    for(int j = frozen_c/2; j < doccpi; ++j)
                    {
                        if(c != d)
                        {
                            if(fabs(epsilon_a[c] - epsilon_a[d]) > 1e-8)
                            {
                                double temp1;
                                double temp2;
                                double temp3;

                                temp1 = 1.0 / (epsilon_a[d] - epsilon_a[c]) * amp_t_dsrg_aa[four_idx(i, j, a, c, nmo)] * amp_t_dsrg_aa[four_idx(i, j, a, d, nmo)] * (epsilon_ijab_aa[four_idx(i, j, a, c, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_aa[four_idx(i, j, a, d, nmo)] * epsilon_ijab_aa[four_idx(i, j, a, d, nmo)])) / (1.0 - pow(e, -S_const * epsilon_ijab_aa[four_idx(i, j, a, c, nmo)] * epsilon_ijab_aa[four_idx(i, j, a, c, nmo)])) - epsilon_ijab_aa[four_idx(i, j, a, d, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_aa[four_idx(i, j, a, c, nmo)] * epsilon_ijab_aa[four_idx(i, j, a, c, nmo)])) / (1.0 - pow(e, -S_const * epsilon_ijab_aa[four_idx(i, j, a, d, nmo)] * epsilon_ijab_aa[four_idx(i, j, a, d, nmo)])));
                                temp2 = 1.0 / (epsilon_a[d] - epsilon_a[c]) * amp_t_dsrg_bb[four_idx(i, j, a, c, nmo)] * amp_t_dsrg_bb[four_idx(i, j, a, d, nmo)] * (epsilon_ijab_bb[four_idx(i, j, a, c, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_bb[four_idx(i, j, a, d, nmo)] * epsilon_ijab_bb[four_idx(i, j, a, d, nmo)])) / (1.0 - pow(e, -S_const * epsilon_ijab_bb[four_idx(i, j, a, c, nmo)] * epsilon_ijab_bb[four_idx(i, j, a, c, nmo)])) - epsilon_ijab_bb[four_idx(i, j, a, d, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_bb[four_idx(i, j, a, c, nmo)] * epsilon_ijab_bb[four_idx(i, j, a, c, nmo)])) / (1.0 - pow(e, -S_const * epsilon_ijab_bb[four_idx(i, j, a, d, nmo)] * epsilon_ijab_bb[four_idx(i, j, a, d, nmo)])));
                                temp3 = 2.0 / (epsilon_a[d] - epsilon_a[c]) * amp_t_dsrg_ab[four_idx(i, j, a, c, nmo)] * amp_t_dsrg_ab[four_idx(i, j, a, d, nmo)] * (epsilon_ijab_ab[four_idx(i, j, a, c, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_ab[four_idx(i, j, a, d, nmo)] * epsilon_ijab_ab[four_idx(i, j, a, d, nmo)])) / (1.0 - pow(e, -S_const * epsilon_ijab_ab[four_idx(i, j, a, c, nmo)] * epsilon_ijab_ab[four_idx(i, j, a, c, nmo)])) - epsilon_ijab_ab[four_idx(i, j, a, d, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_ab[four_idx(i, j, a, c, nmo)] * epsilon_ijab_ab[four_idx(i, j, a, c, nmo)])) / (1.0 - pow(e, -S_const * epsilon_ijab_ab[four_idx(i, j, a, d, nmo)] * epsilon_ijab_ab[four_idx(i, j, a, d, nmo)])));
                                Z_MP2->add(0, 2*d, 2*c, temp1 + temp3);
                                Z_MP2->add(0, 2*d+1, 2*c+1, temp2 + temp3);
                            }
                            else
                            {
                                double temp1;
                                double temp2;
                                double temp3;

                                temp1 = mo_ints_aa[four_idx(i, j, a, c, nmo)] * mo_ints_aa[four_idx(i, j, a, d, nmo)] * (-4.0 * S_const * pow(e, -2.0 * S_const * epsilon_ijab_aa[four_idx(i, j, a, c, nmo)] * epsilon_ijab_aa[four_idx(i, j, a, c, nmo)]) + (1.0 - pow(e, -2.0 * S_const * epsilon_ijab_aa[four_idx(i, j, a, c, nmo)] * epsilon_ijab_aa[four_idx(i, j, a, c, nmo)])) / epsilon_ijab_aa[four_idx(i, j, a, c, nmo)] / epsilon_ijab_aa[four_idx(i, j, a, c, nmo)]);
                                temp2 = mo_ints_bb[four_idx(i, j, a, c, nmo)] * mo_ints_bb[four_idx(i, j, a, d, nmo)] * (-4.0 * S_const * pow(e, -2.0 * S_const * epsilon_ijab_bb[four_idx(i, j, a, c, nmo)] * epsilon_ijab_bb[four_idx(i, j, a, c, nmo)]) + (1.0 - pow(e, -2.0 * S_const * epsilon_ijab_bb[four_idx(i, j, a, c, nmo)] * epsilon_ijab_bb[four_idx(i, j, a, c, nmo)])) / epsilon_ijab_bb[four_idx(i, j, a, c, nmo)] / epsilon_ijab_bb[four_idx(i, j, a, c, nmo)]);
                                temp3 = 2.0 * mo_ints_ab[four_idx(i, j, a, c, nmo)] * mo_ints_ab[four_idx(i, j, a, d, nmo)] * (-4.0 * S_const * pow(e, -2.0 * S_const * epsilon_ijab_ab[four_idx(i, j, a, c, nmo)] * epsilon_ijab_ab[four_idx(i, j, a, c, nmo)]) + (1.0 - pow(e, -2.0 * S_const * epsilon_ijab_ab[four_idx(i, j, a, c, nmo)] * epsilon_ijab_ab[four_idx(i, j, a, c, nmo)])) / epsilon_ijab_ab[four_idx(i, j, a, c, nmo)] / epsilon_ijab_ab[four_idx(i, j, a, c, nmo)]);
                                Z_MP2->add(0, 2*d, 2*c, temp1 + temp3);
                                Z_MP2->add(0, 2*d+1, 2*c+1, temp2 + temp3);
                            }
                        }
                    }
                }
            }
        }
    }  




    // for(int n = frozen_c; n < 2 * doccpi; ++n)
    // {
    //     for(int N = 0; N < frozen_c; ++N)
    //     {
    //         double temp1 = 0.0;

    //         for(int a = 2 * doccpi; a < nso - frozen_v; ++a)
    //         {
    //             for(int b = 2 * doccpi; b < nso - frozen_v; ++b)
    //             {
    //                 for(int j = frozen_c; j < 2 * doccpi; ++j)
    //                 {
    //                     temp1 += so_ints[four_idx(N, j, a, b, nso)] * amp_t_dsrg[four_idx(n, j, a, b, nso)] * (1.0 + pow(e, -S_const * epsilon_ijab[four_idx(n, j, a, b, nso)] * epsilon_ijab[four_idx(n, j, a, b, nso)]));        
    //                 }
    //             }
    //         }
    //         Z_MP2->set(0, n, N, temp1 /(epsilon[n]-epsilon[N]));
    //         Z_MP2->set(0, N, n, Z_MP2->get(0, n, N));
    //     }
    // }
    /********************************** TEST (DONE) **********************************/
    for(int n = frozen_c/2; n < doccpi; ++n)
    {
        for(int N = 0; N < frozen_c/2; ++N)
        {
            double temp1;
            double temp2;
            double temp3;

            for(int a = doccpi; a < nmo - frozen_v/2; ++a)
            {
                for(int b = doccpi; b < nmo - frozen_v/2; ++b)
                {
                    for(int j = frozen_c/2; j < doccpi; ++j)
                    {
                        temp1 = mo_ints_aa[four_idx(N, j, a, b, nmo)] * amp_t_dsrg_aa[four_idx(n, j, a, b, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_aa[four_idx(n, j, a, b, nmo)] * epsilon_ijab_aa[four_idx(n, j, a, b, nmo)]));        
                        temp2 = mo_ints_bb[four_idx(N, j, a, b, nmo)] * amp_t_dsrg_bb[four_idx(n, j, a, b, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_bb[four_idx(n, j, a, b, nmo)] * epsilon_ijab_bb[four_idx(n, j, a, b, nmo)]));        
                        temp3 = 2.0 * mo_ints_ab[four_idx(N, j, a, b, nmo)] * amp_t_dsrg_ab[four_idx(n, j, a, b, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_ab[four_idx(n, j, a, b, nmo)] * epsilon_ijab_ab[four_idx(n, j, a, b, nmo)]));        
                        Z_MP2->add(0, 2*n, 2*N, (temp1 + temp3) /(epsilon_a[n]-epsilon_a[N]));
                        Z_MP2->add(0, 2*n+1, 2*N+1, (temp2 + temp3) /(epsilon_a[n]-epsilon_a[N]));
                    }
                }
            }
            Z_MP2->set(0, 2*N, 2*n, Z_MP2->get(0, 2*n, 2*N));
            Z_MP2->set(0, 2*N+1, 2*n+1, Z_MP2->get(0, 2*n+1, 2*N+1));
        }
    }


    // for(int d = 2 * doccpi; d < nso - frozen_v; ++d)
    // {
    //     for(int D = nso - frozen_v; D < nso; ++D)
    //     {
    //         double temp1 = 0.0;

    //         for(int a = 2 * doccpi; a < nso - frozen_v; ++a)
    //         {
    //             for(int i = frozen_c; i < 2 * doccpi; ++i)
    //             {
    //                 for(int j = frozen_c; j < 2 * doccpi; ++j)
    //                 {
    //                     temp1 += so_ints[four_idx(i, j, a, D, nso)] * amp_t_dsrg[four_idx(i, j, a, d, nso)] * (1.0 + pow(e, -S_const * epsilon_ijab[four_idx(i, j, a, d, nso)] * epsilon_ijab[four_idx(i, j, a, d, nso)]));        
    //                 }
    //             }
    //         }        
    //         Z_MP2->set(0, d, D, temp1/(epsilon[d]-epsilon[D]));
    //         Z_MP2->set(0, D, d, Z_MP2->get(0, d, D));
    //     }
    // } 
    // Z_temp->copy(Z_MP2);
    /********************************** TEST (DONE) **********************************/
    for(int d = doccpi; d < nmo - frozen_v/2; ++d)
    {
        for(int D = nmo - frozen_v/2; D < nmo; ++D)
        {
            double temp1;
            double temp2;
            double temp3;

            for(int a = doccpi; a < nmo - frozen_v/2; ++a)
            {
                for(int i = frozen_c/2; i < doccpi; ++i)
                {
                    for(int j = frozen_c/2; j < doccpi; ++j)
                    {
                        temp1 = mo_ints_aa[four_idx(i, j, a, D, nmo)] * amp_t_dsrg_aa[four_idx(i, j, a, d, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_aa[four_idx(i, j, a, d, nmo)] * epsilon_ijab_aa[four_idx(i, j, a, d, nmo)]));        
                        temp2 = mo_ints_bb[four_idx(i, j, a, D, nmo)] * amp_t_dsrg_bb[four_idx(i, j, a, d, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_bb[four_idx(i, j, a, d, nmo)] * epsilon_ijab_bb[four_idx(i, j, a, d, nmo)]));        
                        temp3 = 2.0 * mo_ints_ab[four_idx(i, j, a, D, nmo)] * amp_t_dsrg_ab[four_idx(i, j, a, d, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_ab[four_idx(i, j, a, d, nmo)] * epsilon_ijab_ab[four_idx(i, j, a, d, nmo)]));        
                        Z_MP2->add(0, 2*d, 2*D, (temp1 + temp3) / (epsilon_a[d] - epsilon_a[D]));
                        Z_MP2->add(0, 2*d+1, 2*D+1, (temp2 + temp3) / (epsilon_a[d] - epsilon_a[D]));
                    }
                }
            }        
            Z_MP2->set(0, 2*D, 2*d, Z_MP2->get(0, 2*d, 2*D));
            Z_MP2->set(0, 2*D+1, 2*d+1, Z_MP2->get(0, 2*d+1, 2*D+1));
        }
    } 




    Z_temp->copy(Z_MP2);







// for( int times = 0; times < 100; ++times)
// {
//     for(int c = 2 * doccpi; c < nso - frozen_v; ++c)
//     {
//         for(int n = frozen_c; n < 2 *doccpi; ++n)
//         {
//             double T_1 = 0.0; 
//             double T_2 = 0.0; 
//             double T_3 = 0.0; 
//             double T_4 = 0.0; 
//             double T_5 = 0.0; 

//             for(int j = frozen_c; j < 2 * doccpi; ++j)
//             {
//                 for(int a = 2 * doccpi; a < nso - frozen_v; ++a)
//                 {
//                     for(int b = 2 * doccpi; b < nso - frozen_v; ++b)
//                     {

//                         T_1 += so_ints[four_idx(c, j, a, b, nso)] * amp_t_dsrg[four_idx(n, j, a, b, nso)] * (1.0 + pow(e, -S_const * epsilon_ijab[four_idx(n, j, a, b, nso)] * epsilon_ijab[four_idx(n, j, a, b, nso)])); 
//                     }
//                 }
//             }

//             for(int i = frozen_c; i < 2 * doccpi; ++i)
//             {
//                 for(int j = frozen_c; j < 2 * doccpi; ++j)
//                 {
//                     for(int a = 2 * doccpi; a < nso - frozen_v; ++a)
//                     {
//                         T_2 -= so_ints[four_idx(i, j, a, n, nso)] * amp_t_dsrg[four_idx(i, j, a, c, nso)] * (1.0 + pow(e, -S_const * epsilon_ijab[four_idx(i, j, a, c, nso)] * epsilon_ijab[four_idx(i, j, a, c, nso)])); 
//                     }
//                 }
//             }

//             for(int i = frozen_c; i < 2 * doccpi; ++i)
//             {
//                 for(int j = frozen_c; j < 2 * doccpi; ++j)
//                 {
//                     for(int a = 2 * doccpi; a < nso - frozen_v; ++a)
//                     {
//                         for(int b = 2 * doccpi; b < nso - frozen_v; ++b)
//                         {
//                             T_4 -= (so_ints[four_idx(i, c, i, n, nso)] - so_ints[four_idx(a, c, a, n, nso)]) * amp_t_dsrg[four_idx(i, j, a, b, nso)] * amp_t_dsrg[four_idx(i, j, a, b, nso)] * (1.0 + pow(e, -S_const * epsilon_ijab[four_idx(i, j, a, b, nso)] * epsilon_ijab[four_idx(i, j, a, b, nso)])) / (1.0 - pow(e, -S_const * epsilon_ijab[four_idx(i, j, a, b, nso)] * epsilon_ijab[four_idx(i, j, a, b, nso)]));
//                         }
//                     }
//                 }
//             }

//             for(int i = frozen_c; i < 2 * doccpi; ++i)
//             {
//                 for(int j = frozen_c; j < 2 * doccpi; ++j)
//                 {
//                     for(int a = 2 * doccpi; a < nso - frozen_v; ++a)
//                     {
//                         for(int b = 2 * doccpi; b < nso - frozen_v; ++b)
//                         {
//                             T_5 += 4.0 * S_const*(so_ints[four_idx(i, c, i, n, nso)] - so_ints[four_idx(a, c, a, n, nso)]) * so_ints[four_idx(i, j, a, b, nso)] * so_ints[four_idx(i, j, a, b, nso)] * pow(e, -2.0 * S_const * epsilon_ijab[four_idx(i, j, a, b, nso)] * epsilon_ijab[four_idx(i, j, a, b, nso)]);
//                         }
//                     }
//                 }
//             }

//             for(int p = 0; p < nso; ++p)
//             {
//                 for(int q = 0; q < nso; ++q)
//                 {
//                     if (p != q)
//                     {
//                         T_3 += so_ints[four_idx(p, n, q, c, nso)] * Z_MP2->get(0, q, p);
//                     }
//                 }
//             }

//             Z_temp->set(0, n, c, (T_1 + T_2 + T_3 + T_4 + T_5) / (epsilon[n] - epsilon[c]));
//             Z_temp->set(0, c, n, Z_temp->get(0, n, c));
//         }
//     }

//     /***********        Z {IA} {AI}         ***********/
//     for(int I = 0; I < frozen_c; ++I)
//     {
//         for(int A = nso - frozen_v; A < nso; ++A)
//         {
//             double temp1 = 0.0;
//             double temp2 = 0.0;
//             double temp3 = 0.0;

//             for(int p = 0; p < nso; ++p)
//             {
//                 for(int q = 0; q < nso; ++q)
//                 {
//                     if (p != q)
//                     {
//                         temp1 += so_ints[four_idx(p, A, q, I, nso)] * Z_MP2->get(0, q, p);
//                     }
//                 }
//             }

//             for(int i = frozen_c; i < 2 * doccpi; ++i)
//             {
//                 for(int j = frozen_c; j < 2 * doccpi; ++j)
//                 {
//                     for(int a = 2 * doccpi; a < nso - frozen_v; ++a)
//                     {
//                         for(int b = 2 * doccpi; b < nso - frozen_v; ++b)
//                         {
//                             temp2 -= ( so_ints[four_idx(i, I, i, A, nso)] - so_ints[four_idx(a, I, a, A, nso)] ) * amp_t_dsrg[four_idx(i, j, a, b, nso)] * amp_t_dsrg[four_idx(i, j, a, b, nso)] * (1.0 + pow(e, -S_const * epsilon_ijab[four_idx(i, j, a, b, nso)] * epsilon_ijab[four_idx(i, j, a, b, nso)])) / (1.0 - pow(e, -S_const * epsilon_ijab[four_idx(i, j, a, b, nso)] * epsilon_ijab[four_idx(i, j, a, b, nso)]));
//                         }
//                     }
//                 }
//             }

//             for(int i = frozen_c; i < 2 * doccpi; ++i)
//             {
//                 for(int j = frozen_c; j < 2 * doccpi; ++j)
//                 {
//                     for(int a = 2 * doccpi; a < nso - frozen_v; ++a)
//                     {
//                         for(int b = 2 * doccpi; b < nso - frozen_v; ++b)
//                         {
//                             temp3 += 4.0 * S_const * so_ints[four_idx(i, j, a, b, nso)] * so_ints[four_idx(i, j, a, b, nso)] * pow(e, -2.0 * S_const * epsilon_ijab[four_idx(i, j, a, b, nso)] * epsilon_ijab[four_idx(i, j, a, b, nso)]) * ( so_ints[four_idx(i, I, i, A, nso)] - so_ints[four_idx(a, I, a, A, nso)] );  
//                         }
//                     }
//                 }
//             }  
//             Z_temp->set(0, I, A, (temp1 + temp2 + temp3) / ( epsilon[I] - epsilon[A] ));
//             Z_temp->set(0, A, I, Z_temp->get(0, I, A));
//         }
//     }        

//     **********        Z {cN} {Nc}         **********
//     for(int c = 2 * doccpi ; c < nso - frozen_v; ++c)
//     {
//         for(int N = 0; N < frozen_c; ++N)
//         {
//             double temp1 = 0.0;
//             double temp2 = 0.0;
//             double temp3 = 0.0;
//             double temp4 = 0.0;

//             for(int p = 0; p < nso; ++p)
//             {
//                 for(int q = 0; q < nso; ++q)
//                 {
//                     if (p != q)
//                     {
//                         temp1 += so_ints[four_idx(p, c, q, N, nso)] * Z_MP2->get(0, q, p);
//                     }
//                 }
//             }

//             for(int i = frozen_c; i < 2 * doccpi; ++i)
//             {
//                 for(int j = frozen_c; j < 2 * doccpi; ++j)
//                 {
//                     for(int b = 2 * doccpi; b < nso - frozen_v; ++b)
//                     {
//                         temp2 -= so_ints[four_idx(i, j, N, b, nso)] * amp_t_dsrg[four_idx(i, j, c, b, nso)] * (1.0 + pow(e, -S_const * epsilon_ijab[four_idx(i, j, b, c, nso)] * epsilon_ijab[four_idx(i, j, b, c, nso)]));
//                     }
//                 }
//             }

//             for(int i = frozen_c; i < 2 * doccpi; ++i)
//             {
//                 for(int j = frozen_c; j < 2 * doccpi; ++j)
//                 {
//                     for(int a = 2 * doccpi; a < nso - frozen_v; ++a)
//                     {
//                         for(int b = 2 * doccpi; b < nso - frozen_v; ++b)
//                         {
//                             temp3 -= ( so_ints[four_idx(i, N, i, c, nso)] - so_ints[four_idx(a, N, a, c, nso)] ) * amp_t_dsrg[four_idx(i, j, a, b, nso)] * amp_t_dsrg[four_idx(i, j, a, b, nso)] * (1.0 + pow(e, -S_const * epsilon_ijab[four_idx(i, j, a, b, nso)] * epsilon_ijab[four_idx(i, j, a, b, nso)])) / (1.0 - pow(e, -S_const * epsilon_ijab[four_idx(i, j, a, b, nso)] * epsilon_ijab[four_idx(i, j, a, b, nso)]));
//                         }
//                     }
//                 }
//             } 

//             for(int i = frozen_c; i < 2 * doccpi; ++i)
//             {
//                 for(int j = frozen_c; j < 2 * doccpi; ++j)
//                 {
//                     for(int a = 2 * doccpi; a < nso - frozen_v; ++a)
//                     {
//                         for(int b = 2 * doccpi; b < nso - frozen_v; ++b)
//                         {
//                             temp4 += 4.0 * S_const * so_ints[four_idx(i, j, a, b, nso)] * so_ints[four_idx(i, j, a, b, nso)] * pow(e, -2.0 * S_const * epsilon_ijab[four_idx(i, j, a, b, nso)] * epsilon_ijab[four_idx(i, j, a, b, nso)]) * ( so_ints[four_idx(i, N, i, c, nso)] - so_ints[four_idx(a, N, a, c, nso)] ); 
//                         }
//                     }
//                 }
//             }
//             Z_temp->set(0, N, c, (temp1 + temp2 + temp3 + temp4) / ( epsilon[N] - epsilon[c] ));
//             Z_temp->set(0, c, N, Z_temp->get(0, N, c));
//         }
//     }   

//     /***********        Z {Cn} {nC}         ***********/
//     for(int C = nso - frozen_v ; C < nso; ++C)
//     {
//         for(int n = frozen_c; n < 2 * doccpi; ++n)
//         {
//             double temp1 = 0.0;
//             double temp2 = 0.0;
//             double temp3 = 0.0;
//             double temp4 = 0.0;

//             for(int p = 0; p < nso; ++p)
//             {
//                 for(int q = 0; q < nso; ++q)
//                 {
//                     if (p != q)
//                     {
//                         temp1 += so_ints[four_idx(p, n, q, C, nso)] * Z_MP2->get(0, q, p);
//                     }
//                 }
//             }

//             for(int j = frozen_c; j < 2 * doccpi; ++j)
//             {
//                 for(int a = 2 * doccpi; a < nso - frozen_v; ++a)
//                 {
//                     for(int b = 2 * doccpi; b < nso - frozen_v; ++b)
//                     {
//                         temp2 += so_ints[four_idx(C, j, a, b, nso)] * amp_t_dsrg[four_idx(n, j, a, b, nso)] * (1.0 + pow(e, -S_const * epsilon_ijab[four_idx(n, j, a, b, nso)] * epsilon_ijab[four_idx(n, j, a, b, nso)]));
//                     }
//                 }
//             }

//             for(int i = frozen_c; i < 2 * doccpi; ++i)
//             {
//                 for(int j = frozen_c; j < 2 * doccpi; ++j)
//                 {
//                     for(int a = 2 * doccpi; a < nso - frozen_v; ++a)
//                     {
//                         for(int b = 2 * doccpi; b < nso - frozen_v; ++b)
//                         {
//                             temp3 += ( so_ints[four_idx(a, n, a, C, nso)] - so_ints[four_idx(i, n, i, C, nso)] ) * amp_t_dsrg[four_idx(i, j, a, b, nso)] * amp_t_dsrg[four_idx(i, j, a, b, nso)] * (1.0 + pow(e, -S_const * epsilon_ijab[four_idx(i, j, a, b, nso)] * epsilon_ijab[four_idx(i, j, a, b, nso)])) / (1.0 - pow(e, -S_const * epsilon_ijab[four_idx(i, j, a, b, nso)] * epsilon_ijab[four_idx(i, j, a, b, nso)])); 
//                         }
//                     }
//                 }
//             }    

//             for(int i = frozen_c; i < 2 * doccpi; ++i)
//             {
//                 for(int j = frozen_c; j < 2 * doccpi; ++j)
//                 {
//                     for(int a = 2 * doccpi; a < nso - frozen_v; ++a)
//                     {
//                         for(int b = 2 * doccpi; b < nso - frozen_v; ++b)
//                         {
//                             temp4 += 4.0 * S_const * so_ints[four_idx(i, j, a, b, nso)] * so_ints[four_idx(i, j, a, b, nso)] * pow(e, -2.0 * S_const * epsilon_ijab[four_idx(i, j, a, b, nso)] * epsilon_ijab[four_idx(i, j, a, b, nso)]) * ( so_ints[four_idx(i, C, i, n, nso)] - so_ints[four_idx(a, C, a, n, nso)] ); 
//                         }
//                     }
//                 }
//             }
//             Z_temp->set(0, n, C, (temp1 + temp2 + temp3 + temp4) / ( epsilon[n] - epsilon[C] ));
//             Z_temp->set(0, C, n, Z_temp->get(0, n, C));
//         }
//     }  
//     Z_MP2->copy(Z_temp);
// }
/********************************** TEST **********************************/



    std::vector<double> Xi_a(nmo, 0.0);
    std::vector<double> Xi_b(nmo, 0.0);
    std::vector<double> Xa_a(nmo, 0.0);
    std::vector<double> Xa_b(nmo, 0.0);
    std::vector<double> Yi_a(nmo, 0.0);
    std::vector<double> Yi_b(nmo, 0.0);
    std::vector<double> Ya_a(nmo, 0.0);
    std::vector<double> Ya_b(nmo, 0.0);

    for(int i = frozen_c/2; i < doccpi; ++i)
    {
        for(int j = frozen_c/2; j < doccpi; ++j)
        {
            for(int a = doccpi; a < nmo - frozen_v/2; ++a)
            {
                for(int b = doccpi; b < nmo - frozen_v/2; ++b)
                {
                    Xi_a[i] += amp_t_dsrg_aa[four_idx(i, j, a, b, nmo)] * amp_t_dsrg_aa[four_idx(i, j, a, b, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_aa[four_idx(i, j, a, b, nmo)] * epsilon_ijab_aa[four_idx(i, j, a, b, nmo)])) / (1.0 - pow(e, -S_const * epsilon_ijab_aa[four_idx(i, j, a, b, nmo)] * epsilon_ijab_aa[four_idx(i, j, a, b, nmo)]));
                    Xi_a[i] += 2.0 * amp_t_dsrg_ab[four_idx(i, j, a, b, nmo)] * amp_t_dsrg_ab[four_idx(i, j, a, b, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)] * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)])) / (1.0 - pow(e, -S_const * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)] * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)]));
                    Xi_b[i] += amp_t_dsrg_bb[four_idx(i, j, a, b, nmo)] * amp_t_dsrg_bb[four_idx(i, j, a, b, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_bb[four_idx(i, j, a, b, nmo)] * epsilon_ijab_bb[four_idx(i, j, a, b, nmo)])) / (1.0 - pow(e, -S_const * epsilon_ijab_bb[four_idx(i, j, a, b, nmo)] * epsilon_ijab_bb[four_idx(i, j, a, b, nmo)]));
                    Xi_b[i] += 2.0 * amp_t_dsrg_ab[four_idx(i, j, a, b, nmo)] * amp_t_dsrg_ab[four_idx(i, j, a, b, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)] * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)])) / (1.0 - pow(e, -S_const * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)] * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)]));
                    Xa_a[a] += amp_t_dsrg_aa[four_idx(i, j, a, b, nmo)] * amp_t_dsrg_aa[four_idx(i, j, a, b, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_aa[four_idx(i, j, a, b, nmo)] * epsilon_ijab_aa[four_idx(i, j, a, b, nmo)])) / (1.0 - pow(e, -S_const * epsilon_ijab_aa[four_idx(i, j, a, b, nmo)] * epsilon_ijab_aa[four_idx(i, j, a, b, nmo)]));
                    Xa_a[a] += 2.0 * amp_t_dsrg_ab[four_idx(i, j, a, b, nmo)] * amp_t_dsrg_ab[four_idx(i, j, a, b, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)] * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)])) / (1.0 - pow(e, -S_const * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)] * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)]));
                    Xa_b[a] += amp_t_dsrg_bb[four_idx(i, j, a, b, nmo)] * amp_t_dsrg_bb[four_idx(i, j, a, b, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_bb[four_idx(i, j, a, b, nmo)] * epsilon_ijab_bb[four_idx(i, j, a, b, nmo)])) / (1.0 - pow(e, -S_const * epsilon_ijab_bb[four_idx(i, j, a, b, nmo)] * epsilon_ijab_bb[four_idx(i, j, a, b, nmo)]));
                    Xa_b[a] += 2.0 * amp_t_dsrg_ab[four_idx(i, j, a, b, nmo)] * amp_t_dsrg_ab[four_idx(i, j, a, b, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)] * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)])) / (1.0 - pow(e, -S_const * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)] * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)]));   
                    
                    Yi_a[i] += mo_ints_aa[four_idx(i, j, a, b, nmo)] * mo_ints_aa[four_idx(i, j, a, b, nmo)] * pow(e, -2.0 * S_const * epsilon_ijab_aa[four_idx(i, j, a, b, nmo)] * epsilon_ijab_aa[four_idx(i, j, a, b, nmo)]);
                    Yi_a[i] += 2.0 * mo_ints_ab[four_idx(i, j, a, b, nmo)] * mo_ints_ab[four_idx(i, j, a, b, nmo)] * pow(e, -2.0 * S_const * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)] * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)]);
                    Yi_b[i] += mo_ints_bb[four_idx(i, j, a, b, nmo)] * mo_ints_bb[four_idx(i, j, a, b, nmo)] * pow(e, -2.0 * S_const * epsilon_ijab_bb[four_idx(i, j, a, b, nmo)] * epsilon_ijab_bb[four_idx(i, j, a, b, nmo)]);
                    Yi_b[i] += 2.0 * mo_ints_ab[four_idx(i, j, a, b, nmo)] * mo_ints_ab[four_idx(i, j, a, b, nmo)] * pow(e, -2.0 * S_const * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)] * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)]);
                    Ya_a[a] += mo_ints_aa[four_idx(i, j, a, b, nmo)] * mo_ints_aa[four_idx(i, j, a, b, nmo)] * pow(e, -2.0 * S_const * epsilon_ijab_aa[four_idx(i, j, a, b, nmo)] * epsilon_ijab_aa[four_idx(i, j, a, b, nmo)]);
                    Ya_a[a] += 2.0 * mo_ints_ab[four_idx(i, j, a, b, nmo)] * mo_ints_ab[four_idx(i, j, a, b, nmo)] * pow(e, -2.0 * S_const * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)] * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)]);
                    Ya_b[a] += mo_ints_bb[four_idx(i, j, a, b, nmo)] * mo_ints_bb[four_idx(i, j, a, b, nmo)] * pow(e, -2.0 * S_const * epsilon_ijab_bb[four_idx(i, j, a, b, nmo)] * epsilon_ijab_bb[four_idx(i, j, a, b, nmo)]);
                    Ya_b[a] += 2.0 * mo_ints_ab[four_idx(i, j, a, b, nmo)] * mo_ints_ab[four_idx(i, j, a, b, nmo)] * pow(e, -2.0 * S_const * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)] * epsilon_ijab_ab[four_idx(i, j, a, b, nmo)]);                   
                }
            }
        }
    }

for( int times = 0; times < 100; ++times)
{
    /***********        Z {nc} {cn} (DONE)         ***********/
    for(int c = doccpi; c < nmo - frozen_v/2; ++c)
    {
        for(int n = frozen_c/2; n < doccpi; ++n)
        {
            double T1_temp1 = 0.0, T1_temp2 = 0.0, T1_temp3 = 0.0;
            double T2_temp1 = 0.0, T2_temp2 = 0.0, T2_temp3 = 0.0;
            double T3_temp1 = 0.0, T3_temp2 = 0.0, T3_temp3 = 0.0;
            double T4_temp1 = 0.0, T4_temp2 = 0.0;
            double T5_temp1 = 0.0, T5_temp2 = 0.0; 

            for(int j = frozen_c/2; j < doccpi; ++j)
            {
                for(int a = doccpi; a < nmo - frozen_v/2; ++a)
                {
                    for(int b = doccpi; b < nmo - frozen_v/2; ++b)
                    {
                        T1_temp1 += mo_ints_aa[four_idx(c, j, a, b, nmo)] * amp_t_dsrg_aa[four_idx(n, j, a, b, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_aa[four_idx(n, j, a, b, nmo)] * epsilon_ijab_aa[four_idx(n, j, a, b, nmo)])); 
                        T1_temp2 += mo_ints_bb[four_idx(c, j, a, b, nmo)] * amp_t_dsrg_bb[four_idx(n, j, a, b, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_bb[four_idx(n, j, a, b, nmo)] * epsilon_ijab_bb[four_idx(n, j, a, b, nmo)])); 
                        T1_temp3 += 2.0 * mo_ints_ab[four_idx(c, j, a, b, nmo)] * amp_t_dsrg_ab[four_idx(n, j, a, b, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_ab[four_idx(n, j, a, b, nmo)] * epsilon_ijab_ab[four_idx(n, j, a, b, nmo)])); 
                    }
                }
            }

            for(int i = frozen_c/2; i < doccpi; ++i)
            {
                for(int j = frozen_c/2; j < doccpi; ++j)
                {
                    for(int a = doccpi; a < nmo - frozen_v/2; ++a)
                    {
                        T2_temp1 -= mo_ints_aa[four_idx(i, j, a, n, nmo)] * amp_t_dsrg_aa[four_idx(i, j, a, c, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_aa[four_idx(i, j, a, c, nmo)] * epsilon_ijab_aa[four_idx(i, j, a, c, nmo)])); 
                        T2_temp2 -= mo_ints_bb[four_idx(i, j, a, n, nmo)] * amp_t_dsrg_bb[four_idx(i, j, a, c, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_bb[four_idx(i, j, a, c, nmo)] * epsilon_ijab_bb[four_idx(i, j, a, c, nmo)])); 
                        T2_temp3 -= 2.0 * mo_ints_ab[four_idx(i, j, a, n, nmo)] * amp_t_dsrg_ab[four_idx(i, j, a, c, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_ab[four_idx(i, j, a, c, nmo)] * epsilon_ijab_ab[four_idx(i, j, a, c, nmo)])); 
                    }
                }
            }


            for(int i = frozen_c/2; i < doccpi; ++i)
            {
                T4_temp1 -= mo_ints_aa[four_idx(i, c, i, n, nmo)] * Xi_a[i];
                T4_temp1 -= mo_ints_ab[four_idx(i, c, i, n, nmo)] * Xi_b[i];
                T4_temp2 -= mo_ints_bb[four_idx(i, c, i, n, nmo)] * Xi_b[i];
                T4_temp2 -= mo_ints_ab[four_idx(i, c, i, n, nmo)] * Xi_a[i];

                T5_temp1 += 4.0 * S_const * mo_ints_aa[four_idx(i, c, i, n, nmo)] * Yi_a[i];
                T5_temp1 += 4.0 * S_const * mo_ints_ab[four_idx(i, c, i, n, nmo)] * Yi_b[i];
                T5_temp2 += 4.0 * S_const * mo_ints_bb[four_idx(i, c, i, n, nmo)] * Yi_b[i];
                T5_temp2 += 4.0 * S_const * mo_ints_ab[four_idx(i, c, i, n, nmo)] * Yi_a[i];
            }

            for(int a = doccpi; a < nmo - frozen_v/2; ++a)
            {
                T4_temp1 += mo_ints_aa[four_idx(a, c, a, n, nmo)] * Xa_a[a];
                T4_temp1 += mo_ints_ab[four_idx(a, c, a, n, nmo)] * Xa_b[a];
                T4_temp2 += mo_ints_bb[four_idx(a, c, a, n, nmo)] * Xa_b[a];
                T4_temp2 += mo_ints_ab[four_idx(a, c, a, n, nmo)] * Xa_a[a];

                T5_temp1 -= 4.0 * S_const * mo_ints_aa[four_idx(a, c, a, n, nmo)] * Ya_a[a];
                T5_temp1 -= 4.0 * S_const * mo_ints_ab[four_idx(a, c, a, n, nmo)] * Ya_b[a];
                T5_temp2 -= 4.0 * S_const * mo_ints_bb[four_idx(a, c, a, n, nmo)] * Ya_b[a];
                T5_temp2 -= 4.0 * S_const * mo_ints_ab[four_idx(a, c, a, n, nmo)] * Ya_a[a];
            }


            for(int p = 0; p < nmo; ++p)
            {
                for(int q = 0; q < nmo; ++q)
                {
                    if (p != q)
                    {
                        T3_temp1 += mo_ints_aa[four_idx(p, n, q, c, nmo)] * Z_MP2->get(0, 2*q, 2*p);
                        T3_temp1 += mo_ints_ab[four_idx(p, n, q, c, nmo)] * Z_MP2->get(0, 2*q+1, 2*p+1);
                        T3_temp2 += mo_ints_bb[four_idx(p, n, q, c, nmo)] * Z_MP2->get(0, 2*q+1, 2*p+1);
                        T3_temp2 += mo_ints_ab[four_idx(p, n, q, c, nmo)] * Z_MP2->get(0, 2*q, 2*p);
                    }
                }
            }

            Z_temp->set(0, 2*n, 2*c, (T1_temp1 + T1_temp3 + T2_temp1 + T2_temp3 + T3_temp1 + T4_temp1 + T5_temp1) / (epsilon_a[n] - epsilon_a[c]));
            Z_temp->set(0, 2*c, 2*n, Z_temp->get(0, 2*n, 2*c));
            Z_temp->set(0, 2*n+1, 2*c+1, (T1_temp2 + T1_temp3 + T2_temp2 + T2_temp3 + T3_temp2 + T4_temp2 + T5_temp2) / (epsilon_a[n] - epsilon_a[c]));
            Z_temp->set(0, 2*c+1, 2*n+1, Z_temp->get(0, 2*n+1, 2*c+1));
        }
    }

    /***********        Z {IA} {AI} (DONE)        ***********/
    for(int I = 0; I < frozen_c/2; ++I)
    {
        for(int A = nmo - frozen_v/2; A < nmo; ++A)
        {
            double T1_temp1 = 0.0, T1_temp2 = 0.0;
            double T2_temp1 = 0.0, T2_temp2 = 0.0;
            double T3_temp1 = 0.0, T3_temp2 = 0.0;

            for(int p = 0; p < nmo; ++p)
            {
                for(int q = 0; q < nmo; ++q)
                {
                    if (p != q)
                    {
                        T1_temp1 += mo_ints_aa[four_idx(p, A, q, I, nmo)] * Z_MP2->get(0, 2*q, 2*p);
                        T1_temp1 += mo_ints_ab[four_idx(p, A, q, I, nmo)] * Z_MP2->get(0, 2*q+1, 2*p+1);
                        T1_temp2 += mo_ints_bb[four_idx(p, A, q, I, nmo)] * Z_MP2->get(0, 2*q+1, 2*p+1);
                        T1_temp2 += mo_ints_ab[four_idx(p, A, q, I, nmo)] * Z_MP2->get(0, 2*q, 2*p);
                    }
                }
            }

            for(int i = frozen_c/2; i < doccpi; ++i)
            {
                T2_temp1 -= mo_ints_aa[four_idx(i, I, i, A, nmo)] * Xi_a[i];
                T2_temp1 -= mo_ints_ab[four_idx(i, I, i, A, nmo)] * Xi_b[i];
                T2_temp2 -= mo_ints_bb[four_idx(i, I, i, A, nmo)] * Xi_b[i];
                T2_temp2 -= mo_ints_ab[four_idx(i, I, i, A, nmo)] * Xi_a[i];

                T3_temp1 += 4.0 * S_const * mo_ints_aa[four_idx(i, I, i, A, nmo)] * Yi_a[i];
                T3_temp1 += 4.0 * S_const * mo_ints_ab[four_idx(i, I, i, A, nmo)] * Yi_b[i];
                T3_temp2 += 4.0 * S_const * mo_ints_bb[four_idx(i, I, i, A, nmo)] * Yi_b[i];
                T3_temp2 += 4.0 * S_const * mo_ints_ab[four_idx(i, I, i, A, nmo)] * Yi_a[i];
            }

            for(int a = doccpi; a < nmo - frozen_v/2; ++a)
            {
                T2_temp1 += mo_ints_aa[four_idx(a, I, a, A, nmo)] * Xa_a[a];
                T2_temp1 += mo_ints_ab[four_idx(a, I, a, A, nmo)] * Xa_b[a];
                T2_temp2 += mo_ints_bb[four_idx(a, I, a, A, nmo)] * Xa_b[a];
                T2_temp2 += mo_ints_ab[four_idx(a, I, a, A, nmo)] * Xa_a[a];

                T3_temp1 -= 4.0 * S_const * mo_ints_aa[four_idx(a, I, a, A, nmo)] * Ya_a[a];
                T3_temp1 -= 4.0 * S_const * mo_ints_ab[four_idx(a, I, a, A, nmo)] * Ya_b[a];
                T3_temp2 -= 4.0 * S_const * mo_ints_bb[four_idx(a, I, a, A, nmo)] * Ya_b[a];
                T3_temp2 -= 4.0 * S_const * mo_ints_ab[four_idx(a, I, a, A, nmo)] * Ya_a[a];
            }

            Z_temp->set(0, 2*I, 2*A, (T1_temp1 + T2_temp1 + T3_temp1) / (epsilon_a[I] - epsilon_a[A]));
            Z_temp->set(0, 2*A, 2*I, Z_temp->get(0, 2*I, 2*A));
            Z_temp->set(0, 2*I+1, 2*A+1, (T1_temp2 + T2_temp2 + T3_temp2) / (epsilon_a[I] - epsilon_a[A]));
            Z_temp->set(0, 2*A+1, 2*I+1, Z_temp->get(0, 2*I+1, 2*A+1));
        }
    }        

    /***********        Z {cN} {Nc} (DONE)        ***********/
    for(int c = doccpi ; c < nmo - frozen_v/2; ++c)
    {
        for(int N = 0; N < frozen_c/2; ++N)
        {
            double T1_temp1 = 0.0, T1_temp2 = 0.0;
            double T2_temp1 = 0.0, T2_temp2 = 0.0;
            double T3_temp1 = 0.0, T3_temp2 = 0.0;
            double T4_temp1 = 0.0, T4_temp2 = 0.0, T4_temp3 = 0.0;

            for(int p = 0; p < nmo; ++p)
            {
                for(int q = 0; q < nmo; ++q)
                {
                    if (p != q)
                    {
                        T1_temp1 += mo_ints_aa[four_idx(p, c, q, N, nmo)] * Z_MP2->get(0, 2*q, 2*p);
                        T1_temp1 += mo_ints_ab[four_idx(p, c, q, N, nmo)] * Z_MP2->get(0, 2*q+1, 2*p+1);
                        T1_temp2 += mo_ints_bb[four_idx(p, c, q, N, nmo)] * Z_MP2->get(0, 2*q+1, 2*p+1);
                        T1_temp2 += mo_ints_ab[four_idx(p, c, q, N, nmo)] * Z_MP2->get(0, 2*q, 2*p);
                    }
                }
            }

            for(int i = frozen_c/2; i < doccpi; ++i)
            {
                T2_temp1 -= mo_ints_aa[four_idx(i, N, i, c, nmo)] * Xi_a[i];
                T2_temp1 -= mo_ints_ab[four_idx(i, N, i, c, nmo)] * Xi_b[i];
                T2_temp2 -= mo_ints_bb[four_idx(i, N, i, c, nmo)] * Xi_b[i];
                T2_temp2 -= mo_ints_ab[four_idx(i, N, i, c, nmo)] * Xi_a[i];

                T3_temp1 += 4.0 * S_const * mo_ints_aa[four_idx(i, N, i, c, nmo)] * Yi_a[i];
                T3_temp1 += 4.0 * S_const * mo_ints_ab[four_idx(i, N, i, c, nmo)] * Yi_b[i];
                T3_temp2 += 4.0 * S_const * mo_ints_bb[four_idx(i, N, i, c, nmo)] * Yi_b[i];
                T3_temp2 += 4.0 * S_const * mo_ints_ab[four_idx(i, N, i, c, nmo)] * Yi_a[i];
            }

            for(int a = doccpi; a < nmo - frozen_v/2; ++a)
            {
                T2_temp1 += mo_ints_aa[four_idx(a, N, a, c, nmo)] * Xa_a[a];
                T2_temp1 += mo_ints_ab[four_idx(a, N, a, c, nmo)] * Xa_b[a];
                T2_temp2 += mo_ints_bb[four_idx(a, N, a, c, nmo)] * Xa_b[a];
                T2_temp2 += mo_ints_ab[four_idx(a, N, a, c, nmo)] * Xa_a[a];

                T3_temp1 -= 4.0 * S_const * mo_ints_aa[four_idx(a, N, a, c, nmo)] * Ya_a[a];
                T3_temp1 -= 4.0 * S_const * mo_ints_ab[four_idx(a, N, a, c, nmo)] * Ya_b[a];
                T3_temp2 -= 4.0 * S_const * mo_ints_bb[four_idx(a, N, a, c, nmo)] * Ya_b[a];
                T3_temp2 -= 4.0 * S_const * mo_ints_ab[four_idx(a, N, a, c, nmo)] * Ya_a[a];
            }

            for(int i = frozen_c/2; i < doccpi; ++i)
            {
                for(int j = frozen_c/2; j < doccpi; ++j)
                {
                    for(int a = doccpi; a < nmo - frozen_v/2; ++a)
                    {
                        T4_temp1 -= mo_ints_aa[four_idx(i, j, a, N, nmo)] * amp_t_dsrg_aa[four_idx(i, j, a, c, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_aa[four_idx(i, j, a, c, nmo)] * epsilon_ijab_aa[four_idx(i, j, a, c, nmo)])); 
                        T4_temp2 -= mo_ints_bb[four_idx(i, j, a, N, nmo)] * amp_t_dsrg_bb[four_idx(i, j, a, c, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_bb[four_idx(i, j, a, c, nmo)] * epsilon_ijab_bb[four_idx(i, j, a, c, nmo)])); 
                        T4_temp3 -= 2.0 * mo_ints_ab[four_idx(i, j, a, N, nmo)] * amp_t_dsrg_ab[four_idx(i, j, a, c, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_ab[four_idx(i, j, a, c, nmo)] * epsilon_ijab_ab[four_idx(i, j, a, c, nmo)])); 
                    }
                }
            }

            Z_temp->set(0, 2*N, 2*c, (T1_temp1 + T2_temp1 + T3_temp1 + T4_temp1 + T4_temp3) / (epsilon_a[N] - epsilon_a[c]));
            Z_temp->set(0, 2*c, 2*N, Z_temp->get(0, 2*N, 2*c));
            Z_temp->set(0, 2*N+1, 2*c+1, (T1_temp2 + T2_temp2 + T3_temp2 + T4_temp2 + T4_temp3) / (epsilon_a[N] - epsilon_a[c]));
            Z_temp->set(0, 2*c+1, 2*N+1, Z_temp->get(0, 2*N+1, 2*c+1));
        }
    }   

    /***********        Z {Cn} {nC}         ***********/
    for(int C = nmo - frozen_v/2; C < nmo; ++C)
    {
        for(int n = frozen_c/2; n < doccpi; ++n)
        {
            double T1_temp1 = 0.0, T1_temp2 = 0.0;
            double T2_temp1 = 0.0, T2_temp2 = 0.0;
            double T3_temp1 = 0.0, T3_temp2 = 0.0;
            double T4_temp1 = 0.0, T4_temp2 = 0.0, T4_temp3 = 0.0;

            for(int p = 0; p < nmo; ++p)
            {
                for(int q = 0; q < nmo; ++q)
                {
                    if (p != q)
                    {
                        T1_temp1 += mo_ints_aa[four_idx(p, n, q, C, nmo)] * Z_MP2->get(0, 2*q, 2*p);
                        T1_temp1 += mo_ints_ab[four_idx(p, n, q, C, nmo)] * Z_MP2->get(0, 2*q+1, 2*p+1);
                        T1_temp2 += mo_ints_bb[four_idx(p, n, q, C, nmo)] * Z_MP2->get(0, 2*q+1, 2*p+1);
                        T1_temp2 += mo_ints_ab[four_idx(p, n, q, C, nmo)] * Z_MP2->get(0, 2*q, 2*p);
                    }
                }
            }

            for(int i = frozen_c/2; i < doccpi; ++i)
            {
                T2_temp1 -= mo_ints_aa[four_idx(i, n, i, C, nmo)] * Xi_a[i];
                T2_temp1 -= mo_ints_ab[four_idx(i, n, i, C, nmo)] * Xi_b[i];
                T2_temp2 -= mo_ints_bb[four_idx(i, n, i, C, nmo)] * Xi_b[i];
                T2_temp2 -= mo_ints_ab[four_idx(i, n, i, C, nmo)] * Xi_a[i];

                T3_temp1 += 4.0 * S_const * mo_ints_aa[four_idx(i, n, i, C, nmo)] * Yi_a[i];
                T3_temp1 += 4.0 * S_const * mo_ints_ab[four_idx(i, n, i, C, nmo)] * Yi_b[i];
                T3_temp2 += 4.0 * S_const * mo_ints_bb[four_idx(i, n, i, C, nmo)] * Yi_b[i];
                T3_temp2 += 4.0 * S_const * mo_ints_ab[four_idx(i, n, i, C, nmo)] * Yi_a[i];
            }

            for(int a = doccpi; a < nmo - frozen_v/2; ++a)
            {
                T2_temp1 += mo_ints_aa[four_idx(a, n, a, C, nmo)] * Xa_a[a];
                T2_temp1 += mo_ints_ab[four_idx(a, n, a, C, nmo)] * Xa_b[a];
                T2_temp2 += mo_ints_bb[four_idx(a, n, a, C, nmo)] * Xa_b[a];
                T2_temp2 += mo_ints_ab[four_idx(a, n, a, C, nmo)] * Xa_a[a];

                T3_temp1 -= 4.0 * S_const * mo_ints_aa[four_idx(a, n, a, C, nmo)] * Ya_a[a];
                T3_temp1 -= 4.0 * S_const * mo_ints_ab[four_idx(a, n, a, C, nmo)] * Ya_b[a];
                T3_temp2 -= 4.0 * S_const * mo_ints_bb[four_idx(a, n, a, C, nmo)] * Ya_b[a];
                T3_temp2 -= 4.0 * S_const * mo_ints_ab[four_idx(a, n, a, C, nmo)] * Ya_a[a];
            }

            for(int j = frozen_c/2; j < doccpi; ++j)
            {
                for(int a = doccpi; a < nmo - frozen_v/2; ++a)
                {
                    for(int b = doccpi; b < nmo - frozen_v/2; ++b)
                    {
                        T4_temp1 += mo_ints_aa[four_idx(C, j, a, b, nmo)] * amp_t_dsrg_aa[four_idx(n, j, a, b, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_aa[four_idx(n, j, a, b, nmo)] * epsilon_ijab_aa[four_idx(n, j, a, b, nmo)])); 
                        T4_temp2 += mo_ints_bb[four_idx(C, j, a, b, nmo)] * amp_t_dsrg_bb[four_idx(n, j, a, b, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_bb[four_idx(n, j, a, b, nmo)] * epsilon_ijab_bb[four_idx(n, j, a, b, nmo)])); 
                        T4_temp3 += 2.0 * mo_ints_ab[four_idx(C, j, a, b, nmo)] * amp_t_dsrg_ab[four_idx(n, j, a, b, nmo)] * (1.0 + pow(e, -S_const * epsilon_ijab_ab[four_idx(n, j, a, b, nmo)] * epsilon_ijab_ab[four_idx(n, j, a, b, nmo)])); 
                    }
                }
            }

            Z_temp->set(0, 2*n, 2*C, (T1_temp1 + T2_temp1 + T3_temp1 + T4_temp1 + T4_temp3) / ( epsilon_a[n] - epsilon_a[C] ));
            Z_temp->set(0, 2*C, 2*n, Z_temp->get(0, 2*n, 2*C));
            Z_temp->set(0, 2*n+1, 2*C+1, (T1_temp2 + T2_temp2 + T3_temp2 + T4_temp2 + T4_temp3) / ( epsilon_a[n] - epsilon_a[C] ));
            Z_temp->set(0, 2*C+1, 2*n+1, Z_temp->get(0, 2*n+1, 2*C+1));
        }
    }  
    Z_MP2->copy(Z_temp);
}






















// Z_MP2->print();

    for(int p = 0; p < nso; ++p)
    {
        for(int q = 0; q < nso; ++q)
        {
            if(p!=q) 
            {
                D_MP2->set(0, p, q, 0.5 * Z_MP2->get(0, p, q));
            }
        }
    }


D_MP2->print();

    double dipole_MP2 = 0.0;

    for(int p = 0; p < nso; ++p)
    {
        for(int q = 0; q < nso; ++q)
        {
            dipole_MP2 += D_MP2->get(0, p, q) * Dp_mo->get(0, p / 2, q / 2);
        }
    }

    // Emp2 = MP2_Energy_SO(eri_mo, F_MO, nso, doccpi, so_ints, epsilon_ijab, frozen_c, frozen_v );
    // Edsrg_pt2 = DSRG_PT2_Energy_SO(eri_mo, F_MO, nso, doccpi, so_ints, epsilon_ijab, S_const, frozen_c, frozen_v );

    /****** test ********/

    Emp2 = MP2_Energy_MO(eri_mo, F_MO, nmo, doccpi, mo_ints_aa, mo_ints_bb, mo_ints_ab, epsilon_ijab_aa, frozen_c, frozen_v);
    Edsrg_pt2 = DSRG_PT2_Energy_MO(eri_mo, F_MO, nmo, doccpi, mo_ints_aa, mo_ints_bb, mo_ints_ab, epsilon_ijab_aa, epsilon_ijab_bb, epsilon_ijab_ab, S_const, frozen_c, frozen_v);


    /****** test ********/    

    //Output
	std::cout << "Energy precision(SCF Iter):   "<< std::setprecision(15) << CVG << std::endl << std::endl;
	// std::cout << "Iteration times:              "<< iternum << std::endl;
    std::cout << "Nuclear repulsion energy:     "<< std::setprecision(15) << Enuc << std::endl;
    std::cout << "Electronic energy:            "<< std::setprecision(15) << Elec << std::endl;
    std::cout << "SCF energy:                   "<< std::setprecision(15) << Escf << std::endl;
    std::cout << "MP2 energy:                   "<< std::setprecision(15) << Emp2 << std::endl;
    std::cout << "Total energy(MP2):            "<< std::setprecision(15) << Escf + Emp2 << std::endl;
    std::cout << "DSRG-PT2 energy:              "<< std::setprecision(15) << Edsrg_pt2 << std::endl;
    std::cout << "Total energy(DSRG-PT2):       "<< std::setprecision(15) << Escf + Edsrg_pt2 << std::endl;
    if(gradient)
    {
        std::cout << "DSRG dipole:                  "<< std::setprecision(15) << dipole_MP2 << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl << std::endl << std::endl;

    for(int i = 0; i < nso; ++i)
    {
        std::cout << i + 1 << "    " << std::setprecision(15) << epsilon[i] << std::endl;
    }
    std::cout << "doccpi   =   " << doccpi << std::endl;
    return ref_wfn;
}
}
} // End namespaces

