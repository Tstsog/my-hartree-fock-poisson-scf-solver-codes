# my-hartree-fock-poisson-scf-solver-codes
## Publisher
Author: Tsogbayar Tsednee (PhD)
Contact: tsog215@gmail.com

## Introduction: 

A matlab code main_he_atom.m computes an eigenenergy for the ground state of helium atom by calling a matlab code he_atom_class.m, which contains a Poisson and self-consistent field calculations and differentiation matrix of a pseudospectral method.  
HF_Poisson_solver_H2_ground_state.m and DFT_Poisson_solver_H2_ground_state.m codes calculate the total energy of ground state of molecular hydrogen, respectively. Both use a Poisson solver to calculate interelectronic interaction.  


## Requirement: 
Any version of Matlab 

## Implementation details and running

The code main_he_atom.m uses the he_atom_class.m. A class he_atom_class calculates the Poisson equation, of which solution is used to compute electron-electron interaction, and Hartree-Fock self-consistent equation and the Legendre differentiation matrix with collocation at the Legendre-Gauss-Lobatto nodes and corresponding weights.    
The codes HF_Poisson_solver_H2_ground_state.m and DFT_Poisson_solver_H2_ground_state.m use pseudospectral method, as well. 

You may download them and run the main_he_atom.m, HF_Poisson_solver_H2_ground_state.m and DFT_Poisson_solver_H2_ground_state.m directly. 

## Copyright / License 

These codes are published as freeware. Basically, you have the right to use them and modify any of them. 

Find the GNU General Public License at:
https://www.gnu.org/licenses/gpl-3.0.en.html
