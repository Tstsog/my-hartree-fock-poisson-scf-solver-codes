%%% Computes eigen energy for the lowest four doubly exited states of helium atom
% Written by Tsogbayar Tsednee (PhD), California State University Northridge;
% Contact: tsog215@gmail.com
% Reference: Ts.Tsogbayar & D. L. Yeager,  Chinese Physics B 26, 083101 (2017)
% March 18, 2020
% 
clear; clc;
Z = 2.;   % charge of Helium nuclie; 
N = 64;   % number of grid point; you may change it
rb = 0.0; % beginning point of coordinate r 
re = 15.; % ending point of coordnate r; you may change it
% unit used here is atomic unit
%%%
itermax = 100; tol = 10^(-12);
%
[r, psi, Vee, En1, En_total] = he_atom_class.he_1s2_state(Z,N,rb,re,itermax,tol);
%
plot(r, r.*psi, 'b');
[Vee, En1, En_total ] % Vee is an interelectronic interaction; 
                      % En1 is 1s orbital energy
                      % En_total is total energy for the 1s2 state

% [Vee, En1, En_total ] = 1.025768869900022  -0.917955562856104  -2.861679995612231