function [Q] = c_NR(param, fns)
% Newton-Raphson scheme to solve the dEL equations
%
% :param param: initial conditions, boundary conditions, material parameter
% :param fns: symbolic functions of residual and tangent
%
% :returns: discrete trajectory of the beam

%---------------------------
%           1      2      3       ...  i
%    Q = q_0   q_1   q_2   ...  q_n
%---------------------------
Q = zeros(param.dim_q,length(param.time)); % initial Q to store trajectory
Q(:,1) = param.q_ref;

%% ++++1. step:  sympletic momentum==========================================
%        known: q_0,  solve: q_1
%---------------------------

% index of free dof, i.e. dof without boundary conditions
% in 1. step, fix mechanical dof, solve maxwell equations only
id = [];
switch param.eboun
    case 'end'
        id = 10:param.dim_q - 6*param.n_kn_nodes-3; % free dof
    otherwise
        inn = 0;
        for in = 2:param.n_kn_nodes
            if(rem(in-1,param.nec) == 0)
            else
                id(inn+1:inn+3) = param.iu.bm(in).ov; % free dof
                inn = inn+3;
            end
        end
end

du = zeros(param.dim_q - 6*param.n_kn_nodes, 1); % initialize the generalized configuration eq.(78)
q_1 = Q(:,1); % initialize q_1 with q_0
q_1 = update(du, q_1, param); % nodal reparametrization eq.(79), impose b.c.

% residual and tangent
[F, dF] = RK0(q_1, param.q_ref, param, fns);

% null space projection
T_d = null_space_matrix(param.q_ref, param); % null space matrix P(q_n) in eq.(80)
F = T_d'*F; % residual 
T_D = null_space_matrix(q_1, param); % da/du in eq.(80)
Kt = T_d'* dF *T_D; % tangent

% crossing out rows and colums where b.c. are impose
Kt_new = Kt(id, id);
F_new = F(id);
%---------------------NR loop to solve q1------------------------------
iter = 0;
while norm(F_new) > param.NR.nr_tol
    if iter > param.NR.max_iter
        error('too many NR iterations')
    end
    disp(['step-0 ' '  iter-',num2str(iter), '  Resi-',num2str(norm(F_new))])
     
    %------------------------------
    du_new = - Kt_new\F_new;
    du = zeros(param.dim_q - 6*param.n_kn_nodes, 1);
    du(id) = du_new;
    %------------------------------
    
    q_1 = update(du, q_1, param);
    
    % Resisum and Tangent
    [F, dF] = RK0(q_1, param.q_ref, param, fns);
    
    F = T_d'*F;
    T_D = null_space_matrix(q_1, param);
    Kt = T_d'* dF *T_D;
    
    Kt_new = Kt(id, id);
    F_new = F(id);
    
    iter = iter +1;
end
Q(:,2) = q_1;
%-----------------------NR loop end----------------------------------

%% +++++++ from 2. steps:  variational integrator===================================
%       known: q_n-1, q_n    solve: q_n+1       (n=1,2,3...)
%-------------------------------------------------

% index of free dof, i.e. dof without boundary conditions
id = [];
switch param.eboun
    case 'end'
        id = 10:param.dim_q - 6*param.n_kn_nodes-3; % no boundary, free dof
    otherwise
        inn = 0;
        for in = 2:param.n_kn_nodes
            if(rem(in-1,param.nec) == 0)
                id(inn+1:inn+6) = [param.iu.bm(in).wphi, param.iu.bm(in).om]; % mech. dof free
                inn = inn+6;
            else
                id(inn+1:inn+9) = [param.iu.bm(in).wphi, param.iu.bm(in).om, param.iu.bm(in).ov]; % all free dof
                inn = inn+9;
            end
        end
end

for i = 2:length(param.time)-1
    q_nm1 = Q(:, i-1);
    q_n = Q(:, i);
    q_np1 = q_n; % initial q_n+1 with q_n
    
    % residual and tangent
    [F, dF] = RK(q_np1, q_nm1, q_n, param, fns);
    
    % null space projection
    T_d = null_space_matrix(q_n, param); % null space matrix P(q_n) in eq.(80)
    F = T_d'*F; % residual
    T_D = null_space_matrix(q_np1, param); % da/du in eq.(80)
    Kt = T_d'* dF *T_D; % tangent
    
    % crossing out rows and colums where b.c. are impose
    Kt_new = Kt(id, id);
    F_new = F(id);
    
    %---------------------NR loop to solve q_n+1------------------------------
    iter = 0;
    while norm(F_new) > param.NR.nr_tol
        if iter > param.NR.max_iter
            error('too many NR iterations')
        end
        
        disp(['step-', num2str(i), '  iter-',num2str(iter), '  Resi-', num2str(norm(F_new))])
        
        %------------------------------
        du_new = - Kt_new\F_new;
        du = zeros(param.dim_q - 6*param.n_kn_nodes, 1);
        du(id) = du_new;
        %------------------------------
        
        q_np1 = update(du, q_np1, param); % nodal reparametrization eq.(79)
        
        [F, dF] = RK(q_np1, q_nm1, q_n, param, fns);
        
        F = T_d'*F;
        T_D = null_space_matrix(q_np1, param);
        Kt = T_d'* dF *T_D;
        
        Kt_new = Kt(id, id);
        F_new = F(id);
        
        iter = iter +1;
    end
    Q(:,i+1) = q_np1;
    %-----------------------NR loop end----------------------------------
end

end
