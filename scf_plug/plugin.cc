/*
 * @BEGIN LICENSE
 *
 * scf_plug by Psi4 Developer, a plugin to:
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
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libpsio/psio.hpp"
#include <math.h>
#include <iomanip>
#include <vector>

namespace psi{ namespace scf_plug {

extern "C"
int read_options(std::string name, Options& options)
{
    if (name == "SCF_PLUG"|| options.read_globals()) {
        options.add_int("PRINT", 1);
        options.add_double("CVG", 1);
    }

    return true;
}

void FormDensityMatrix(SharedMatrix D, SharedMatrix C, int dim, int occ){

    D->zero();

    for(int i=0; i<dim; ++i){

		for(int j=0; j<dim; ++j){

			for(int k=0; k<occ; ++k){

				D->add(0, i, j, C->get(0, i, k) * C->get(0, j, k));

			}
		}
	}
}

void FormNewFockMatrix(SharedMatrix F, SharedMatrix H, SharedMatrix D, SharedMatrix eri, int dim){

    for(int p=0; p<dim; ++p){

    	for(int q=0; q<dim; ++q){

    		F->set(0, p, q, H->get(0, p, q));

    		for(int r=0; r<dim; ++r){

    			for(int s=0; s<dim; ++s){

    				F->add(0, p, q, D->get(0, r, s) * (2.0*eri->get(0, p*dim+q, r*dim+s) - eri->get(0, p*dim+r, q*dim+s)));

    			}
    		}
    	}
    }
}

double ElecEnergy(double Elec, SharedMatrix D, SharedMatrix H, SharedMatrix F, int dim){

    Elec = 0;

    for(int p=0; p<dim; ++p){

    	for(int q=0; q<dim; ++q){

    		Elec+= D->get(0, p, q)*(H->get(0, p, q)+F->get(0, p, q));

    	}
    }

    return Elec;
}

extern "C"
SharedWavefunction scf_plug(SharedWavefunction ref_wfn, Options& options)
{

//Parameters Declaration
    std::shared_ptr<Molecule> molecule = Process::environment.molecule();
    molecule->update_geometry();
    molecule->print();

    std::shared_ptr<BasisSet> ao_basisset=ref_wfn->basisset();
    MintsHelper mints(ref_wfn);

    int print = options.get_int("PRINT");
    double cvg = options.get_double("CVG");
    int dims[]={ao_basisset->nbf()};
    double CVG = options.get_double("CVG");
    int iternum = 1;
    double energy_pre;
    int doccpi = 0;
    double Enuc = molecule->nuclear_repulsion_energy();
    double Elec, Etot;
    int irrep_num = ref_wfn->nirrep();

    std::shared_ptr<MatrixFactory> factory(new MatrixFactory);
    factory->init_with(1,dims,dims);
  
    SharedMatrix overlap = mints.ao_overlap();
    SharedMatrix kinetic = mints.ao_kinetic();
    SharedMatrix potential = mints.ao_potential();
	SharedMatrix Omega = factory->create_shared_matrix("Omega");
    SharedMatrix F (new Matrix("Fock matrix", 1, dims, dims, 0));
    SharedMatrix C (new Matrix("C matrix", 1, dims, dims, 0));
    SharedMatrix D (new Matrix("Density matrix", 1, dims, dims, 0));
    SharedMatrix S (new Matrix("S matrix", 1, dims, dims, 0));
    SharedMatrix H = factory->create_shared_matrix("H");
    SharedMatrix eri = mints.ao_eri();
    SharedMatrix evecs (new Matrix("evecs", 1, dims, dims, 0));
    SharedVector evals (new Vector("evals", 1, dims));


    Dimension doccpi_add = ref_wfn->doccpi();

    for(int i=0;i<irrep_num;i++){
    	doccpi += doccpi_add[i];
    }

    std::cout<<std::endl;

//Create H matrix
    H->copy(kinetic);
    H->add(potential);

//Create S^(-1/2) Matrix

    Omega->zero();
    overlap->diagonalize(evecs, evals);

    for(int i=0; i < dims[0]; ++i){
    	Omega->set(0, i, i, 1.0/sqrt(evals->get(0,i)));
    }
     
    S = Matrix::triplet(evecs, Omega, evecs, false, false, true);

//Create original fock matrix using transformation on H

	F = Matrix::triplet(S, H, S, true, false, false);

//Create C matrix

    F->diagonalize(evecs, evals);
    C = Matrix::doublet(S, evecs, false, false);

//Create Density matrix

    FormDensityMatrix(D, C, dims[0], doccpi);

//Create new Fock matrix

	FormNewFockMatrix(F, H, D, eri, dims[0]);

//Calculate the energy

	Elec = ElecEnergy(Elec, D, H, F, dims[0]);
    Etot = Elec + Enuc;
    energy_pre = Etot; 

    std::cout<<"Nuclear repulsion energy:     "<< std::setprecision(15)<< Enuc<< std::endl;
    std::cout<<"Electronic energy:            "<< std::setprecision(15)<< Elec<< std::endl;
    std::cout<<"Total energy:                 "<< std::setprecision(15)<< Etot<< std::endl<< std::endl;

    Etot += CVG + 1.0;

//SCF iteration

while( fabs(energy_pre - Etot) > CVG){

	energy_pre = Etot;

    F = Matrix::triplet(S, F, S, true, false, false);
    F->diagonalize(evecs, evals);
    C = Matrix::doublet(S, evecs, false, false);

    FormDensityMatrix(D, C, dims[0], doccpi);
	FormNewFockMatrix(F, H, D, eri, dims[0]);

    Elec = ElecEnergy(Elec, D, H, F, dims[0]);
    Etot = Elec + Enuc;

    iternum++;
}
//Output
	std::cout<<"Energy precision:             "<< std::setprecision(15)<< CVG<< std::endl<< std::endl;
	std::cout<<"Iteration times:              "<< iternum<< std::endl;
    std::cout<<"Nuclear repulsion energy:     "<< std::setprecision(15)<< Enuc<< std::endl;
    std::cout<<"Electronic energy:            "<< std::setprecision(15)<< Elec<< std::endl;
    std::cout<<"Total energy:                 "<< std::setprecision(15)<< Etot<< std::endl<< std::endl;

    return ref_wfn;
}
}} // End namespaces

