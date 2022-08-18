function [W_e] = int_potEn_cm(eref, eq, param)
% strain energy computed from kinematics in contunuum mechanics
%
% :param eref: nodal configuraitons of beam element in reference configuraiton
% :param eq: nodal configuraitons of beam element in current configuraiton
% :param param: parameters
%
% :returns: strain energy

%% nodal values
phi0_node =[eref(1:3,1),eref(1:3,2)]; % 3*2  phi 1,2
phi_node   =[  eq(1:3,1),  eq(1:3,2)];

dir0_node =[eref(4:12,1),eref(4:12,2)]; % 9*2   d
dir_node   =[  eq(4:12,1),  eq(4:12,2)];

v0_node   =[eref(13:15,1),eref(13:15,2)]; % 3*2, v-electric potential
v_node     =[  eq(13:15,1),  eq(13:15,2)];

% Gauss integration
ngp=2; % 2 Gauss points are used
g1=1/sqrt(3); % 1/sqrt(3), see PW book table 4.1
w1=1;
xsi=[-g1; g1]; wp=[w1; w1];

%% ------------ integration over arc length --------------
W_mech=0; W_ele=0;

for i=1:ngp
    
    % Lagrange type linear shape functions:, see PW book eq(4.17)
    N(1,1)=(1-xsi(i))/2; %N1
    N(1,2)=(1+xsi(i))/2;%N2       1*2
    
    % dN/dxi
    dNr(1,1)=-1/2; %dN1dxi    1*2
    dNr(1,2)= 1/2;  %dN2dxi
    
    % d_phi/d_xi
    dp0_dx=phi0_node(:,1)'*dNr(1) + phi0_node(:,2)'*dNr(2); % Jacobian matrix=dX/dxi=XI*dN/dxi, see PW book eq(4.35)
    dp_dx =phi_node(:,1)'*dNr(1) + phi_node(:,2)'*dNr(2); % 1*3
    
    % d(xi)
    dir0=N*dir0_node'; % 1*9
    dir =N*dir_node';
    
    % dd/dxi
    ddir0_dx=dNr*dir0_node'; % 1*9
    ddir_dx =dNr*dir_node';
    
    % v - electric potential
    vv = N*v_node'; %1*3
    
    % dv/dxi  - electric
    dv0_dx=dNr*v0_node'; % 1*3
    dv_dx =dNr*v_node';
    
    p1=[1:3]; p2=[4:6]; p3=[7:9];
        
    % ds/dxi = ||d_phi0/dxi || : Jacobian
    ndp_dx=norm(dp0_dx);
    
    if ndp_dx<10*eps
        disp('length equal or less than zero!')
    end
    
    % d_phi/d_s = d_phi/d_xi * dxi/ds
    dp0_ds  =  dp0_dx/ndp_dx; % 1*3
    dp_ds    =   dp_dx/ndp_dx;
    
    % d_dir/d_s = d_dir/d_s * dxi/ds
    ddir0_ds = ddir0_dx/ndp_dx; % 1*9
    ddir_ds   = ddir_dx/ndp_dx;
    
    % dv/ds = d_v/d_xi * dxi/ds
    dv0_ds  = dv0_dx/ndp_dx; % 1*3
    dv_ds    = dv_dx/ndp_dx;
    
    %v-strain
    v=(dp_ds - dir(p3)); % 1*3
    
    % Lmanbda
    Lam = dir(p1)' * dir0(p1) + dir(p2)' * dir0(p2) + dir(p3)' * dir0(p3); % 3*3
    
    %% ----------in tegration over cross section-------------------------------
    W1=0; W2=0;
    
    % gauss integration for last term
    X_node=[-param.b/2  -param.b/2;
                    param.b/2  -param.b/2;
                    param.b/2   param.b/2;
                  -param.b/2   param.b/2]; % nodes of the cross section 4*2
    alpha=[-1/sqrt(3)  1/sqrt(3)  1/sqrt(3) -1/sqrt(3)];
    beta =[-1/sqrt(3) -1/sqrt(3) 1/sqrt(3) 1/sqrt(3)];
    wht=1;
    for pi=1:4
        % shape function
        Na(1,1)=1/4*(1-alpha(pi))*(1-beta(pi));
        Na(1,2)=1/4*(1+alpha(pi))*(1-beta(pi));
        Na(1,3)=1/4*(1+alpha(pi))*(1+beta(pi));
        Na(1,4)=1/4*(1-alpha(pi))*(1+beta(pi));
        XY=Na*X_node;%(x,y)
        
        % derivative of shape function
        dNa=[-1/4*(1- beta(pi))   1/4*(1-beta(pi))    1/4*(1+beta(pi))  -1/4*(1+beta(pi));...
                 -1/4*(1-alpha(pi)) -1/4*(1+alpha(pi))  1/4*(1+alpha(pi))  1/4*(1-alpha(pi))];
        
        Jac=0;
        for nn=1:4
            Jac = Jac + X_node(nn,:)'*dNa(:,nn)'; %dX/d_eta
        end
        Je=det(Jac);
        
        % electricla field
        Ef = - ( vv(2)*dir0(p1) + vv(3)*dir0(p2) + (dv_ds(1) + XY(1)*dv_ds(2) + XY(2)*dv_ds(3) )*dir0(p3) ); % 1*3
        
        % deformation gradient
        a = dp_ds - dir(p3) + XY(1)*ddir_ds(p1) + XY(2)*ddir_ds(p2);
        F =  (eye(3,3) + a' * dir(p3)) * Lam;
        
        J=1+dot(a,dir(p3));
        
        A_inv = eye(3,3) - (a' * dir(p3))/J;
        
        F_inv = Lam' * A_inv;
        
        % right Cauchy Green tensor
        C = F' * F;
        C_inv = F_inv * F_inv';
        
        % mechanical free energy,  integration over cross section
        W1 = W1 + (param.mu/2*(trace(C)-3) - param.mu*log(J) +param.lam/2*log(J)*log(J) )*Je*wht;
        %W1 = W1 + (param.mu/2*(trace(C)-3) - param.mu*dot(a,dir(p3)) +param.lam/2*dot(a,dir(p3))*dot(a,dir(p3)) )*Je*wht;
        
        % double contraction of tensors
        CEE=0; CinvEE=0;
        EdE=Ef' * Ef;
        for ii=1:3
            for jj=1:3
                CEE = CEE + C(ii,jj)*EdE(ii,jj);
                CinvEE = CinvEE + C_inv(ii,jj)*EdE(ii,jj);
            end
        end
        
        % electrical free energy,  integration over cross section
        W2 = W2 + (param.c1*dot(Ef,Ef) + param.c2*CEE - param.ep0/2*J*CinvEE )*Je*wht;
        
    end
    %% ---------------------------------------------
    
    W_mech = W_mech+ W1*ndp_dx*wp(i); % integration over center line
    
    W_ele = W_ele+ W2*ndp_dx*wp(i); % integration over center line
    
end

W_e = W_mech + W_ele; % total free energy
end

