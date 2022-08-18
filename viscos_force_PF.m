function [P_vis, F_e, ndp_dx] = viscos_force_PF(eref, eq, e_dot, XY, param)
% 1st Piola-Kirchhoff stress (Kelvin–Voigt viscosity model) and deformation gradient F in element j
%
% :param eref: nodal configuraitons of beam element in reference configuraiton
% :param eq: nodal configuraitons of beam element in current configuraiton
% :param e_dot: nodal velocity of beam element in current configuraiton
% :param XY: coordinates of integration point in cross section
% :param param: parameters
%
% :returns: 1st Piola-Kirchhoff stress and deformation gradient F

% nodal values
phi0_node      =[ eref(1:3,1),   eref(1:3,2)]; % 3*2
phi_node        =[   eq(1:3,1),   eq(1:3,2)];
phi_dot_node =[e_dot(1:3,1),  e_dot(1:3,2)];

dir0_node       =[  eref(4:12,1),   eref(4:12,2)]; % 9*2
dir_node         =[     eq(4:12,1),  eq(4:12,2)];
dir_dot_node   =[e_dot(4:12,1),  e_dot(4:12,2)];

xsi=0.0; % one gauss point is applied!!

% shape functions: N
N(:,1)=(1-xsi)/2;   N(:,2)=(1+xsi)/2; % 1*2

% dN/dxi
dNr(:,1)=-1/2;      dNr(:,2)= 1/2;  %1*2

% d_phi/d_xi  d_phi_dot/d_xi
dp0_dx=dNr*phi0_node'; % 1*3
dp_dx =dNr*phi_node';
dp_dot_dx =dNr*phi_dot_node';

% d(xi)  d_dot(xi)
dir0=N*dir0_node'; % 1*9
dir =N*dir_node';
dir_dot =N*dir_dot_node';

% dd/dxi   dd_dot/dxi
ddir0_dx=dNr*dir0_node'; % 1*9
ddir_dx =dNr*dir_node';
ddir_dot_dx =dNr*dir_dot_node';

p1=[1:3]; p2=[4:6]; p3=[7:9];

%%

% Ja:  ds/dxi = ||d_phi0/dxi || : Jacobian
ndp_dx=norm(dp0_dx);

if ndp_dx<10*eps
    disp('length equal or less than zero!')
end

% d_phi/d_s = d_phi/d_xi * dxi/ds
dp_ds    =   dp_dx/ndp_dx; % 1*3
dp_dot_ds    =   dp_dot_dx/ndp_dx;

% d_dir/d_s = d_dir/d_s * dxi/ds
ddir_ds   = ddir_dx/ndp_dx; % 1*9
ddir_dot_ds   = ddir_dot_dx/ndp_dx;

% Lmanbda, Lmanbda_dot
Lam = dir(p1)' * dir0(p1) + dir(p2)' * dir0(p2) + dir(p3)' * dir0(p3); % 3*3
Lam_dot = dir_dot(p1)' * dir0(p1) + dir_dot(p2)' * dir0(p2) + dir_dot(p3)' * dir0(p3); % 3*3

% evaluated at point: xi, alpha, beta
a = dp_ds - dir(p3) + XY(1)*ddir_ds(p1) + XY(2)*ddir_ds(p2);

a_dot = dp_dot_ds - dir_dot(p3) + XY(1)*ddir_dot_ds(p1) + XY(2)*ddir_dot_ds(p2);

F_e =  (eye(3,3) + a' * dir(p3)) * Lam;

%F_dot = (eye(3,3) + a_dot' * dir(p3)) * Lam_dot; % wrong !!
F_dot = Lam_dot + a_dot'*dir0(p3);

J=1-dot(a,dir(p3));

A_inv = eye(3,3) - (a' * dir(p3))/J;

F_inv = Lam' * A_inv;

C_inv = F_inv * F_inv';

P_vis = J*param.eta*(F_inv' * F_dot' * F_inv' + F_dot * C_inv)/2;
end
