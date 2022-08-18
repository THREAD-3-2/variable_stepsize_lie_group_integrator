function fns = b_dEL_AD(param)
% derive discrete Euler-Lagrange equations and tangents using automatic differentiation
%
% :param param: initial conditions, boundary conditions, material parameter
%
% :returns: symbolic functions of residual and tangent

% call casadi library
import casadi.*

%% Configuration and Lagrange
dim_q = param.dim_q;

% symbolic variables for autodiff
q_n_var = MX.sym('q_n_var', dim_q, 1); % q_n, it is global vector
q_np1_var = MX.sym('q_np1_var', dim_q, 1); % q_n+1
q_nm1_var = MX.sym('q_nm1_var', dim_q, 1); % q_n-1

% discrete Lagrangian
Lminus = LagrangeEqu(q_n_var, q_np1_var, param.q_ref, param); % n, n+1
Lplus  = LagrangeEqu(q_nm1_var, q_n_var, param.q_ref, param); % n-1, n

% dL/dqn
dLminus = jacobian(Lminus, q_n_var)';  %dL(qn, qn+1)/dqn  % direvative to q_n !!!
dLplus = jacobian(Lplus, q_n_var)'; %dL(qn-1, qn)/dqn

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% Viscoelastic forces
f_extminus  =  -1/2*viscos_force(q_n_var, q_np1_var, param.q_ref, param);
f_extplus = -1/2*viscos_force(q_nm1_var, q_n_var, param.q_ref, param);

%% Constraints
g = constraints(q_n_var,param.q_ref, param); % symbolic fucntion
fns.g = casadi.Function('g_fnc', {q_n_var}, {g}); % casadi function, symbolic variables can be imposed with numerical values 

G = jacobian(g,q_n_var); % automatic differentiation by casadi
fns.G = casadi.Function('G_fnc', {q_n_var}, {G}); % save as casadi fucntion

%% momentum
% p_minus_n
pminus_symb = -dLminus - f_extminus;
fns.pminus = casadi.Function('pminus_symb',{q_n_var, q_np1_var}, {pminus_symb});

%% +++++++++++++Variational integrator+++++++++++++++++++++++++++++++++++++++++++++++++++++
% Residual
R_vi_symb = dLplus + dLminus + f_extminus + f_extplus; % dEL in eq.(73)
fns.R_vi = casadi.Function('R_vi', {q_nm1_var, q_n_var, q_np1_var}, {R_vi_symb});

% Tangent matrix
Kq_vi_symb = jacobian(R_vi_symb, q_np1_var);  % dR/dq in eq.(80)
fns.Kq_vi = casadi.Function('Kq_vi', {q_nm1_var, q_n_var, q_np1_var}, {Kq_vi_symb});

%% +++++++++++ sympletic momentum scheme, only used in 1. step ++++++++++++++++++++
% Residual 
R_sm_symb = param.p_0 + dLminus + f_extminus; % Legendre transformation eq.(87)
fns.R_sm = casadi.Function('R_sm', {q_n_var, q_np1_var}, {R_sm_symb});

% Tangent matrix
Kq_sm_symb = jacobian(R_sm_symb, q_np1_var); % dR/dq in eq.(80)
fns.Kq_sm = casadi.Function('Kq_sm', {q_n_var, q_np1_var}, {Kq_sm_symb});

end