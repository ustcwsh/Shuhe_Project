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

#include <vector>


namespace psi{ namespace scf_plug {

extern "C"
int read_options(std::string name, Options& options)
{
    if (name == "SCF_PLUG"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
        options.add_double("CVG", 1);
    }

    return true;
}

extern "C"
SharedWavefunction scf_plug(SharedWavefunction ref_wfn, Options& options)
{
    int print = options.get_int("PRINT");
    double cvg = options.get_double("CVG");

    std::shared_ptr<Molecule> molecule = Process::environment.molecule();
    molecule->update_geometry();
    molecule->print();

    double Enuc = molecule->nuclear_repulsion_energy();
    outfile->Printf("\tNuuclear repulsion energy: %f\n", Enuc);

    std::shared_ptr<BasisSet> ao_basisset=ref_wfn->basisset();

    int dims[]={ao_basisset->nbf()};

    std::shared_ptr<MatrixFactory> factory(new MatrixFactory);
    factory->init_with(1,dims,dims);

    MintsHelper mints(ref_wfn);
    

    SharedMatrix overlap = mints.ao_overlap();
    SharedMatrix kinetic = mints.ao_kinetic();
    SharedMatrix potential = mints.ao_potential();

    SharedMatrix H = factory->create_shared_matrix("H");
    H->copy(kinetic);
    H->add(potential);
    H->print();

    SharedMatrix eri = mints.ao_eri();

    SharedMatrix Omega = factory->create_shared_matrix("Omega");
    Omega->zero();

    SharedVector evals (new Vector("evals", 1, dims));
    SharedMatrix evecs (new Matrix("evecs", 1, dims, dims, 0));
    overlap->diagonalize(evecs,evals);

    for(int i=0; i < dims[0]; i++){
    	Omega->set(0,i,i,1.0/sqrt(evals->get(0,i)));
    }
     

    SharedMatrix S=Matrix::triplet(evecs,Omega,evecs,false,false,true);
    S->print();
   
   





    // Typically you would build a new wavefunction and populate it with data
    return ref_wfn;
}

}} // End namespaces

