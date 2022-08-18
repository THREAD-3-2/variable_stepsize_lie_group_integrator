% main function
%
%  An electromechanically coupled beam, July,2022, Dengpeng Huang @ LTD-FAU
%  Ref.: D. Huang, S. Leyendecker. Computational Mechanics, 69(2022)(3):805-824.


clear all; clc; close all;
tic

% add paths
addpath('casadi-osx-matlabR2015a-v3.4.5'); addpath('pre_post_processing'); addpath('integrator');

% initialization
a_ini;

% Residual R (discrete Euler-Langerage equaitons) and Tangent K
fns = b_dEL_AD(param);

% solve dEL with Newton-Raphson scheme
res.Q = c_NR(param, fns); 
toc

% postprocessing
Hamilton; % plot the Hamiltonia total energy 
Energy % plot the energy of the beam
plot_FE(res.Q(:,end), param, res, 'v') % plot displacement u/ electric potential v at state
plot_u % plot the displacement of the end node
% plot_FE_mov(res.Q, param, el_col)  % movie
