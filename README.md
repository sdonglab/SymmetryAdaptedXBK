# Symmetry-Adapted XBK Method for Electronic Structure Theory on Quantum Annealers

## Methods

This code implements the algorithm detailed in

- Joseph Desroches, Sijia S. Dong*. Electronic Structure Theory with Molecular Point Group Symmetries on Quantum Annealers. arXiv:2502.00235 https://doi.org/10.48550/arXiv.2502.00235

This code calculates the ground-state energies of electronic Hamiltonians using a combination of the Xia-Bian-Kais (XBK) transformation (Xia, R., Bian, T., & Kais, S. (2017). Electronic Structure calculations and the Ising Hamiltonian. The Journal of Physical Chemistry B, 122(13), 3384–3395. https://doi.org/10.1021/acs.jpcb.7b10371) and symmetry-adapted encodings based on the full Boolean symmetry group (Picozzi, D., & Tennyson, J. (2023). Symmetry-adapted encodings for qubit number reduction by point-group and other Boolean symmetries. Quantum Science and Technology, 8(3), 035026. https://doi.org/10.1088/2058-9565/acd86c). 

We provide a new implementation of the XBK method in `utils.py`, based on exact diagonalization. The `XBK_method.py`, `helper_funcs.py`, and `fieldmath.py` files were from Copenhaver's [repo](https://github.com/jcopenh/Quantum-Chemistry-with-Annealers) for the original implementation of the XBK method to perform the mapping from a qubit Hamiltonian to an Ising Hamiltonian. This original implementation accompanies the author's paper (Copenhaver, J., Wasserman, A., & Wehefritz-Kaufmann, B. (2021). Using quantum annealers to calculate ground state properties of molecules. The Journal of Chemical Physics, 154(3). https://doi.org/10.1063/5.0030397).

## Usage

We suggest using a virtual environment; for example, using [Anaconda](https://www.anaconda.com/). We provide a `.yml` file containing the required packages and the versions tested to run the code.

We provide a tutorial `tutorial.ipynb` to guide new users to perform an example calculation using our method and accompanying functions.

##
Copyright (c) 2025, Dong Lab @ Northeastern University
