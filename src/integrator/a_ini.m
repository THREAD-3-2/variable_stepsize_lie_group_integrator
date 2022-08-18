% set Initial configuration, material parameters and boundary conditions
%
%% parameter setting
% boundary conditions
param.eboun = 'end'; % 'end' (contraction) / 'torsion' / 'shear' / 'bending'/'test'

switch param.eboun
    case 'end'
        param.timestep  = 1e-4; % in [ms]
        param.totaltime = 0.1; % in [ms]
        param.L = 0.1;
        param.b = 0.02; param.h = 0.02;
        param.ncell = 1;
        param.nec = 5; % number ele in each DEA cell
        param.n_el_nodes = param.ncell*param.nec; % number of element
        param.n_kn_nodes = param.n_el_nodes+1;% number of node
        param.eta = 0; % viscosity
    case 'torsion'
        param.timestep  = 1e-3; % in [ms]
        param.totaltime = 0.2; % in [ms]
        param.L = 0.1;
        param.b = 0.005; param.h = 0.005;
        param.ncell = 40;
        param.nec = 1; % number ele in each DEA cell
        param.n_el_nodes = param.ncell*param.nec; % number of element
        param.n_kn_nodes = param.n_el_nodes+1;% number of node
        param.eta = 0.9; % viscosity
    case 'shear'
        param.timestep  = 1e-3; % in [ms]
        param.totaltime = 0.5; % in [ms]
        param.L = 0.05;
        param.b = 0.005; param.h = 0.005;
        param.ncell = 20;
        param.nec = 2; % number ele in each DEA cell
        param.n_el_nodes = param.ncell*param.nec; % number of element
        param.n_kn_nodes = param.n_el_nodes+1;% number of node
        param.eta=3; % viscosity
    case 'bending'
        param.timestep  = 1e-3; % in [ms]
        param.totaltime = 0.5; % in [ms]
        param.L = 0.1;
        param.b = 0.005; param.h = 0.005;
        param.ncell = 40;
        param.nec = 1; % number ele in each DEA cell
        param.n_el_nodes = param.ncell*param.nec; % number of element
        param.n_kn_nodes = param.n_el_nodes+1;% number of node
        param.eta = 10; % viscosity
end
param.time = 0:param.timestep:param.totaltime;
param.N_timesteps = length(param.time);
param.dim_q = 15*param.n_kn_nodes; % number of dof

% material parameters, table 1
param.rho = 1; % [g/mm^3]
param.mu = 233; % Mpa Lames constant
param.lam = 999.8;% MPa Lames constant
param.ep0 = 8.854e-12; % C/(Vm)
param.c1 = 5e-8; % N/V^2
param.c2 = 1e-9;  % N/V^2

% gravity
param.g = 9.81e-3; % grav const in [mm/s^2]
param.g_dir = [0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0]; % x y z d1 d2 d3 *phi*

%% initial configutation q
q_ref = path_s(param.L,param);
q_ref_tmp = reshape(q_ref,15,param.n_kn_nodes);
param.length = sqrt(q_ref_tmp(1,end).^2+q_ref_tmp(2,end).^2+q_ref_tmp(3,end).^2);

param.width = param.b*ones(param.n_kn_nodes,1);
param.height = param.h*ones(param.n_kn_nodes,1);

param.A = param.b*param.h;
I1 = param.b^3*param.h/12;
I2 = param.b*param.h^3/12;
param.inertia = [param.A*param.rho,  I1*param.rho,  I2*param.rho];

q_ref = reshape(q_ref,15*param.n_kn_nodes,1);

param.q_ref = q_ref; % reference configuration for strain calculation
param.p_0 = zeros(size(param.q_ref));

%% elemental mass matrix
param.Me = zeros(15*2,15*2,param.n_el_nodes);
param.Me_bar = zeros(9*2,9*2,param.n_el_nodes);
param.Me_bar_inv = zeros(9*2,9*2,param.n_el_nodes);

for ie=1:param.n_el_nodes % *global assembly
    e_ref_i = [q_ref_tmp(:,ie) q_ref_tmp(:,ie+1)]; % two nodes q1_0, q2_0
    param.Me(:,:,ie) = elmass01(e_ref_i,param.inertia);
    
    % reduced mass matrix, see Leyendecker ZAMM eq.114
    param.Me_bar(:,:,ie) = [param.Me(1:9,1:9,ie) param.Me(1:9,16:24,ie);
                                       param.Me(16:24,1:9,ie) param.Me(16:24,16:24,ie)];
    
    param.Me_bar_inv(:,:,ie) = inv(param.Me_bar(:,:,ie));
end

% solver settings
param.NR.max_iter = 100;
param.NR.nr_tol = 1e-10; % newton rapson tolerance

%% beam idices
for ii=1:param.n_kn_nodes
    param.iq.bm(ii).phi = [15*ii-14:15*ii-12];
    param.iq.bm(ii).d1 = [15*ii-11:15*ii-9];
    param.iq.bm(ii).d2 = [15*ii-8:15*ii-6];
    param.iq.bm(ii).d3 = [15*ii-5:15*ii-3];
    param.iq.bm(ii).v = [15*ii-2:15*ii];
    
    % indices of reduced variable u
    param.iu.bm(ii).wphi = [9*ii-8:9*ii-6]; % phi
    param.iu.bm(ii).om = [9*ii-5:9*ii-3]; % omega
    param.iu.bm(ii).ov = [9*ii-2:9*ii]; % v
end
