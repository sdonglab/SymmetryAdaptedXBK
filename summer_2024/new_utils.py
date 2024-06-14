# ------------------------------------------------------------------------------------------------------------------------------------------
# imports
# ------------------------------------------------------------------------------------------------------------------------------------------

import sys
import math
import quantumsymmetry
import pyscf
import openfermion
import numpy as np
import neal

helper_dir = "/home/jwdesroches/python/Ga2QuAMES/symmetry/helper_funcs/"
sys.path.append(helper_dir)

from XBK_method import *
from QCC_method import *
from fieldmath import *

# ------------------------------------------------------------------------------------------------------------------------------------------
# function definitions
# ------------------------------------------------------------------------------------------------------------------------------------------

def get_unsimplified_hamiltonian(atom, basis, charge, spin, method = "RHF", verbose = False):
    """ Placeholder Text"""

    mol = pyscf.gto.Mole()
    mol.atom = atom
    mol.symmetry = False
    mol.basis = basis
    mol.charge = charge
    mol.spin = spin
    mol.verbose = 0
    mol.build()

    if method == "RHF":
        mf = pyscf.scf.RHF(mol)
    elif method == "UHF":
        mf = pyscf.scf.UHF(mol)
    else:
        print(f"Error: {method} method has not been added yet.")
        
    mf.kernel()

    qubit_op, fermion_op = quantumsymmetry.core.get_hamiltonian(mol, mf)
    
    if verbose == True:
        print(f"qubit op: {qubit_op}")
        print(f"\n fermion op: {fermion_op}")
        
    return qubit_op

# ------------------------------------------------------------------------------------------------------------------------------------------

def get_simplified_hamiltonian(atom, basis, charge, spin):
    """ Placeholder Text"""
    
    qubit_op = quantumsymmetry.reduced_hamiltonian(atom, basis, charge, spin, verbose = False)
    return qubit_op

# ------------------------------------------------------------------------------------------------------------------------------------------

def get_ising_info(qubit_op, r):
    """ Placeholder Text"""
    
    m = openfermion.utils.count_qubits(qubit_op)
    
    qubit_Hs, qubit_Cs = [],[]
    for p in range(int(math.ceil(r/2+1))):
        qubit_Hs += [XBK_transform(qubit_op, r, p)]
        qubit_Cs += [construct_C(m, r, p)]

    return qubit_Hs, qubit_Cs

# ------------------------------------------------------------------------------------------------------------------------------------------

def get_energy_annealing(hamiltonian, r=1, starting_lam = 2, num_samples = 1000, strength = 1e3, 
                         real_annealing=False, direct_diagonalization=True):
    """ Placeholder Text"""

    if real_annealing:
        pass # add functionaltiy later
    else:
        sampler =  neal.SimulatedAnnealingSampler() 
    
    if direct_diagonalization:
        if np.shape(openfermion.linalg.qubit_operator_sparse(hamiltonian))!= (2,2):
            true_energy = openfermion.linalg.get_ground_state(openfermion.linalg.qubit_operator_sparse(hamiltonian))
        else:
            true_energy,_ = np.linalg.eig(openfermion.linalg.qubit_operator_sparse(hamiltonian).toarray())
            true_energy = true_energy.real

    else:
        true_energy = None
    
    ising_H, ising_C = get_ising_info(hamiltonian, r)
    annealed_energy, _ = XBK(ising_H, ising_C, r, sampler, starting_lam=0, num_samples=1000, strength=1e3, verbose=False)
    
    return true_energy, annealed_energy
    
# ------------------------------------------------------------------------------------------------------------------------------------------
     
        
    
    