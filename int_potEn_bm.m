function [W_e] = int_potEn_bm(eref, eq, param)
% strain energy computed from beam kinematics
%
% :param eref: nodal configuraitons of beam element in reference configuraiton
% :param eq: nodal configuraitons of beam element in current configuraiton
% :param param: parameters
%
% :returns: strain energy

%% nodal values of beam element
phi0_node = [eref(1:3,1),eref(1:3,2)]; % 3*2  phi
phi_node   = [  eq(1:3,1),  eq(1:3,2)];

dir0_node = [eref(4:12,1),eref(4:12,2)]; % 9*2   d
dir_node   = [  eq(4:12,1),  eq(4:12,2)];

v0_node   = [eref(13:15,1),eref(13:15,2)]; % 3*2, v-electric potential
v_node     = [  eq(13:15,1),  eq(13:15,2)];

%% gauss point
ngp = 2; % number of gauss point in arc length
g1 = 1/sqrt(3); % 1/sqrt(3), see PW book table 4.1
w1 = 1;
xsi = [-g1; g1]; wp = [w1; w1];

%% ------------ integration over arc length --------------
W_mech = 0; W_elec = 0;
for i = 1:ngp
     % shape functions: N, linear, see PW book eq(4.17)
    N(1,1) = (1 - xsi(i))/2; %N1
    N(1,2) = (1 + xsi(i))/2;%N2       1*2
    
    % dN/dxi
    dNr(1,1) = -1/2; %dN1dxi    1*2
    dNr(1,2) = 1/2;  %dN2dxi
    
    % d_phi/d_xi
    dp0_dx = phi0_node(:,1)'*dNr(1) + phi0_node(:,2)'*dNr(2); % Jacobian matrix=dX/dxi=XI*dN/dxi, see PW book eq(4.35)
    dp_dx = phi_node(:,1)'*dNr(1) + phi_node(:,2)'*dNr(2); % 1*3
    
    % d(xi)
    dir0 = N*dir0_node'; % 1*9
    dir  = N*dir_node';
    
    % dd/dxi
    ddir0_dx = dNr*dir0_node'; % 1*9
    ddir_dx = dNr*dir_node';
    
    % v - electric potential
    vv = N*v_node'; %1*3
    
    % dv/dxi  - electric
    dv0_dxi = dNr*v0_node'; % 1*3
    dv_dxi = dNr*v_node';
    
    p1 = [1:3]; p2 = [4:6]; p3 = [7:9]; % index of directors
    
    %------------------------------------
    % numerical integration of energy over beam element
    
    % ds/dxi = ||d_phi0/dxi || : Jacobian
    ndp_dx = norm(dp0_dx); % scalar
    
    % dN/ds
    dN_ds = dNr/ndp_dx; % 1*2
    
    % d_phi/d_s = d_phi/d_xi * dxi/ds
    dphi0_ds  =  dp0_dx/ndp_dx; % 1*3
    dphi_ds    =   dp_dx/ndp_dx;
    
    % d_dir/d_s = d_dir/d_s * dxi/ds
    ddir0_ds = ddir0_dx/ndp_dx; % 1*9
    ddir_ds   = ddir_dx/ndp_dx;
    
    % dv/ds = d_v/d_xi * dxi/ds
    dv0_ds  = dv0_dxi/ndp_dx; % 1*3
    dv_ds    = dv_dxi/ndp_dx;
    
    % Lmanbda, rotation matrix
    Lam = dir(p1)' * dir0(p1) + dir(p2)' * dir0(p2) + dir(p3)' * dir0(p3); % 3*3
    
    % electricla field
    %Ef = - ( vv(2)*dir0(p1) + vv(3)*dir0(p2) + (dv_ds(1) + XY(1)*dv_ds(2) + XY(2)*dv_ds(3) )*dir0(p3)); % 1*3
    %a = dp_ds - dir(p3) + XY(1)*ddir_ds(p1) + XY(2)*ddir_ds(p2);
    
    % beam strain for shear and elongation
    g = [dphi_ds*dir(p1)'-dphi0_ds*dir0(p1)';
        dphi_ds*dir(p2)'-dphi0_ds*dir0(p2)';
        dphi_ds*dir(p3)'-dphi0_ds*dir0(p3)']; % or gama = dp_ds - dir(p3);
    
    g0 = Lam'*g;
    
    % beam strain for bending and torsion
    k = [ddir_ds(p2)*dir(p3)' - ddir_ds(p3)*dir(p2)' - (ddir0_ds(p2)*dir0(p3)' - ddir0_ds(p3)*dir0(p2)');
        ddir_ds(p3)*dir(p1)' - ddir_ds(p1)*dir(p3)' - (ddir0_ds(p3)*dir0(p1)' - ddir0_ds(p1)*dir0(p3)');
        ddir_ds(p1)*dir(p2)' - ddir_ds(p2)*dir(p1)' - (ddir0_ds(p1)*dir0(p2)' - ddir0_ds(p2)*dir0(p1)')]/2;
    
    k0 = Lam'*k;
    
    % geometric quantities
    d1 = dir(p1);
    d2 = dir(p2);
    d3 = dir(p3);
    d01 = dir0(p1);
    d02 = dir0(p2);
    d03 = dir0(p3);
    
    % electric potential and strain-like variable
    p = [-vv(2), -vv(3), -dv_ds(1)]; % p = [-alpha, -beta, -dp/ds]
    q = [dv_ds(2); dv_ds(3); 0]; % q = [dalpha/ds, dbeta/ds, 0]
    
    % strain energy for beam
    % W beam strain gama,kappa
    [W1, W2, W3, W4] = int_potEn_bm_w(param.b, param.h, param.mu, param.lam, param.c1, param.c2, p, q, d1, d2, d3, d01, d02, d03, g, g0, k, k0);
    
    %------------------------------------
    
    W_mech = W_mech+ (W1 + W2 + W3)*ndp_dx*wp(i); % integration over center line
    
    W_elec = W_elec+ W4*ndp_dx*wp(i); % integration over center line
    
end
W_e = W_mech + W_elec; % total free energy
end

