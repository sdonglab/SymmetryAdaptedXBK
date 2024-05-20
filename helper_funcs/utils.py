from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import BravyiKitaevMapper
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.mappers import ParityMapper
from utils import *
import quantumsymmetry
import numpy as np
import matplotlib.pyplot as plt
from pyscf import gto, scf
from quantumsymmetry.core import get_hamiltonian 
from XBK_method import *
import math

def electronic_to_bk_hamiltonian(atom = "H 0 0 0; H 0 0 0.735", basis = "sto3g", charge=0, spin=0, unit=DistanceUnit.ANGSTROM, 
                                 mapper = "bk", display_output=False, display_type=False):
    
    driver = PySCFDriver(atom=atom, basis=basis, charge=charge, spin=spin, unit=unit,)
    fermionic_op = driver.run().hamiltonian.second_q_op()
    if mapper == "bk":        
        qubit_hamiltonian = BravyiKitaevMapper().map(fermionic_op)
        mapper_title = "Bravyi-Kitaev"
    elif mapper == "jw":
        qubit_hamiltonian = JordanWignerMapper().map(fermionic_op)
        mapper_title = "Jordan-Wigner"
    elif mapper == "p":
        qubit_hamiltonian = ParityMapper().map(fermionic_op)
        mapper_title = "Parity"
    
    if display_output is True:
        print(f"The {mapper_title} qubit Hamiltonian for this system is given by: \n {qubit_hamiltonian}")
    if display_type is True:
        print(f"The {mapper_title} qubit Hamitlonian is of type: \n {type(qubit_hamiltonian)}.")
        
    return qubit_hamiltonian

def get_unsimplified_hamiltonian(atom, basis, charge, spin):

    mol = gto.Mole()
    mol.atom = atom
    mol.symmetry = False
    mol.basis = basis
    mol.charge = charge
    mol.spin = spin
    mol.verbose = 0
    mol.build()
    
    mf = scf.RHF(mol)
    mf.kernel()

    qubit_op, fermion_op = get_hamiltonian(mol, mf)
    return qubit_op

def convert_to_ising_hamiltonian(qubit_op, r):
    ising_hamiltonian = []
    
    for p in range(int(math.ceil(r/2+1))):
        ising_hamiltonian += [XBK_transform_prime(qubit_op, r, p)]

    return ising_hamiltonian[0]
