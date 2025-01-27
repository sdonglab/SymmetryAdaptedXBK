# ------------------------------------------------------------------------------------------------------------------------------------------
# imports
# ------------------------------------------------------------------------------------------------------------------------------------------

import math
import quantumsymmetry
import pyscf
import openfermion
import numpy as np
import neal
import scipy
import time

from XBK_method import *
from fieldmath import *
from helper_functions import *

# ------------------------------------------------------------------------------------------------------------------------------------------
# function definitions
# ------------------------------------------------------------------------------------------------------------------------------------------

def get_unsimplified_hamiltonian(atom, basis, charge, spin, method = "RHF", verbose = False):
    """
    This function generates an unsimplified Hamiltonian for a given molecule using the PySCF and quantumsymmetry libraries.

    Parameters:
    atom (str): The atomic structure of the molecule.
    basis (str): The basis set used for the molecular orbitals.
    charge (int): The total charge of the molecule.
    spin (int): The total spin of the molecule.
    method (str, optional): The method used to calculate the molecular orbitals. Defaults to "RHF".
    verbose (bool, optional): If True, prints the qubit and fermion operators. Defaults to False.

    Returns:
    qubit_op: A QubitOperator object representing the unsimplified Hamiltonian.
    """
    
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
    """
    This function generates a simplified Hamiltonian for a given molecule using the quantumsymmetry library.

    Parameters:
    atom (str): The atomic structure of the molecule.
    basis (str): The basis set used for the molecular orbitals.
    charge (int): The total charge of the molecule.
    spin (int): The total spin of the molecule.

    Returns:
    qubit_op: A QubitOperator object representing the simplified Hamiltonian.
    
    Note: 
    This function was created simply to match our convention for get_unsimplified_hamiltonian().
    """
    
    qubit_op = quantumsymmetry.reduced_hamiltonian(atom, basis, charge, spin, verbose = False)
    return qubit_op

# ------------------------------------------------------------------------------------------------------------------------------------------

def calculate_ising_ground_state(qubit_hamiltonian, r_value):
    """
    Calculates the ground state properties of an qubit Hamiltonian using our exact diagonalization implementation of the XBK method. 

    Parameters:
    qubit_hamiltonian (QubitOperator): The qubit Hamiltonian representing the molecule.
    r_value (float): The r parameter for the XBK transformation.

    Returns:
    tuple: A tuple containing the ground state energy, the ground state wavefunction, and the Ising Hamiltonian.
    """
    
    N = openfermion.utils.count_qubits(qubit_hamiltonian)

    ising_energies = []
    ising_ground_states = []
    if r_value == 1:
        ising_H = XBK_transform(qubit_hamiltonian, 1, 1)
        true_ising_hamiltonian = ising_H

        if r_value*N == 1:
            true_ising_energy, true_ground_state = np.linalg.eigh(openfermion.linalg.get_sparse_operator(ising_H).toarray())
            true_ising_energy = true_ising_energy[0]
            true_ground_state = true_ground_state[:,0]
        else: 
            true_ising_energy, true_ground_state = openfermion.linalg.get_ground_state(openfermion.linalg.get_sparse_operator(ising_H))
            
    else:
        ising_energies = []
        ising_ground_states = []

        ising_Hs, ising_Cs = [],[]
        for p in range(int(math.ceil(r_value/2+1))):
            ising_H = XBK_transform(qubit_hamiltonian, r_value, p)
            ising_C = construct_C(N, r_value, p)
            ising_Hs += [ising_H]
            ising_Cs += [ising_C]

            ising_C_matrix = openfermion.linalg.get_sparse_operator(ising_C)

            if r_value*N == 1:

                ising_energy, ising_gs = np.linalg.eigh(openfermion.linalg.get_sparse_operator(ising_H).toarray())
                ising_energy = ising_energy[0]

            else:

                sparse_operator = openfermion.linalg.get_sparse_operator(ising_H)
                ising_energy, ising_gs = openfermion.linalg.get_ground_state(sparse_operator)

                sumBsq = round(np.linalg.norm(ising_C_matrix @ ising_gs))

                if sumBsq == 0:
                    pass
                ising_energy = ising_energy/sumBsq
                ising_energies.append(ising_energy)
                ising_ground_states.append(ising_gs)

        ising_energies = np.array(ising_energies)
        index = int(np.argmin(ising_energies))

        true_ising_hamiltonian = ising_Hs[index]
        true_ising_energy = ising_energies[index]
        true_ground_state = ising_ground_states[index]

    return true_ising_energy, true_ground_state, true_ising_hamiltonian
  
    
# ------------------------------------------------------------------------------------------------------------------------------------------
def return_geometry(mol_name, bond_length):
    """
    This function generates a pySCF geometry string for a given molecule and bond length.

    Parameters:
    mol_name (str): The name of the molecule. The following molecules are supported: "H2", "He2", "LiH", "H2O", "BH3", "NH3", "CH4", "CO", "F2", "Li2", "O2", "N2".
    bond_length (float): The bond length for the molecule.

    Returns:
    str: A geometry string representing the molecule. If the molecule is not supported, returns None.
    """
    
    if mol_name == "H2":
        geom = f"H 0 0 0; H 0 0 {bond_length}"

    elif mol_name == "He2":
        geom = f"He 0 0 0; He 0 0 {bond_length}"

    elif mol_name == "LiH":
        geom = f"Li 0 0 0; H 0 0 {bond_length}"

    elif mol_name == "H2O":
        bond_angle = np.deg2rad(104.5 / 2)  
        geom = (f"O 0 0 0; "
                f"H 0 {round(bond_length * np.sin(bond_angle), 4)} {round(bond_length * np.cos(bond_angle), 4)}; "
                f"H 0 {-round(bond_length * np.sin(bond_angle), 4)} {round(bond_length * np.cos(bond_angle), 4)}")

    elif mol_name == "BH3":
        bond_angle = np.deg2rad(120)  
        geom = (f"B 0 0 0; "
                f"H {bond_length} 0 0; "
                f"H {-round(bond_length * np.cos(bond_angle), 4)} {round(bond_length * np.sin(bond_angle), 4)} 0; "
                f"H {-round(bond_length * np.cos(bond_angle), 4)} {-round(bond_length * np.sin(bond_angle), 4)} 0")

    elif mol_name == "NH3":
        bond_angle = np.deg2rad(107.8)  
        r = round(bond_length * np.sin(bond_angle / 2),4)
        z = round(bond_length * np.cos(bond_angle / 2),4)
        geom = (f"N 0 0 0; "
                f"H {r} {r} {z}; "
                f"H {r} {-r} {z}; "
                f"H {-r} 0 {z}")

    elif mol_name == "CH4":
        geom = (f"C 0 0 0; "
                f"H {bond_length} {bond_length} {bond_length}; "
                f"H {bond_length} {-bond_length} {-bond_length}; "
                f"H {-bond_length} {bond_length} {-bond_length}; "
                f"H {-bond_length} {-bond_length} {bond_length}")
        
    elif mol_name == "CH2":
        bond_angle = np.deg2rad(22.2485793248)
        geom = (f"C 0 0 0; "
                f"H 0 {np.round(bond_length*np.cos(bond_angle),4)} {-np.round(bond_length*np.sin(bond_angle),4)}; "
                f"H 0 {-np.round(bond_length*np.cos(bond_angle))} {np.round(bond_length*np.sin(bond_angle),4)}")

    elif mol_name == "CO":
        geom = f"C 0 0 0; O 0 0 {bond_length}"

    elif mol_name == "F2":
        geom = f"F 0 0 0; F 0 0 {bond_length}"

    elif mol_name == "Li2":
        geom = f"Li 0 0 0; Li 0 0 {bond_length}"

    elif mol_name == "O2":
        geom = f"O 0 0 0; O 0 0 {bond_length}"

    elif mol_name == "N2":
        geom = f"N 0 0 0; N 0 0 {bond_length}"

    else:
        print("This geometry has not been added yet. Please try again.")
        geom = None
    
    return geom

# ------------------------------------------------------------------------------------------------------------------------------------------
# end of file
# ------------------------------------------------------------------------------------------------------------------------------------------