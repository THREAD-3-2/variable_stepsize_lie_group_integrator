function [Me]=elmass01(eref, inertia)
% (consistent) mass matrix for beam element
%
% :param eref: eref = [..., phi^A,(d1)^A,(d2)^A,(d3)^A, ... ](0); A=1,2
% :param inertia: inertia = (rho*A, M_1, M_2)
%
% :returns: Me, consistent elemental mass matrix (15*2 x 15*2)

%-----------------------------------------------------------------------------
elnode=2;

I_3=eye(3);

%% shape functions and its direvatives

phi0_node=eref(1:3, :); % coordinates of nodes, phi

% Gauss integration
ngp=2;
g1=1/sqrt(3); % 1/sqrt(3), see PW book table 4.1
w1=1;
xsi=[-g1; g1]; wp=[w1; w1];

% shape functions: linear, see PW book eq(4.17)
N(:,1)=(1-xsi)/2; N(:,2)=(1+xsi)/2;

% derivatives of shape functions with respect to normalized coordinate
dNr([1:2],1)=-1/2; dNr([1:2],2)= 1/2;

%% mass matrix, see P.Wriggers FEM book p126
% d_phi/d_xi
dp0_dx=dNr*phi0_node'; % Jacobian matrix=dX/dxi=XI*dN/dxi, see PW book eq(4.35)

% elemental mass matrix size =(elnode*eldof)*(elnode*eldof)=30*30
eldof=15*elnode;

% Me = M_11 M_12
%         M_21 M_22       M_AB=m_AB*I_(12*12), A=1,2; B=1,2
Me=zeros(eldof,eldof);

rho_A=inertia(1); M_1=inertia(2); M_2=inertia(3);

% mass matrix:
for i=1:ngp % I-gauss point
    
    ndp_dx=norm(dp0_dx(i,:)); %det(J)=ds/dxi
    
    if ndp_dx<10*eps
        disp('length equal or less than zero!')
    end
    
    mass_kl=N(i,:)'*N(i,:);  %N_A(xi_I)*N_B(xi_I)
    
    for k=1:elnode % A-loop node
        ka=(k-1)*15;
        k_phi=[ka+1: ka+3];
        k_1=  [ka+4: ka+6];
        k_2=  [ka+7: ka+9];
        
        for l=1:elnode % B-loop node
            la=(l-1)*15;
            l_phi=[la+1: la+3];
            l_1=  [la+4: la+6];
            l_2=  [la+7: la+9];
            
            % see Betsch 2006 eq(18), eq(23) and Leye. ZAMM 2008 eq(114)
            
            % M_AB = m_AB*I_(12*12) = [INT*rhoA*I_3     0            0        0];
            %                                        [      0       INT*E1*I_3      0         0];
            %                                        [      0                    INT*E2*I_3    0];
            %                                        [      0          0            0              0];
            
            % mass matrix comes from kenetic energy T, since T is independet
            % of d3, thus last diagnal elenment in mass matrix is 0
            % INT = sum_I N_A(xi_I)*N_B(xi_I)*wp_I*det(J)
            Me(k_phi,l_phi) = Me(k_phi,l_phi) +  rho_A*mass_kl(k,l)*I_3*ndp_dx*wp(i);
            Me(k_1  ,l_1)    = Me(k_1  ,l_1)     +  M_1  *mass_kl(k,l)*I_3*ndp_dx*wp(i);
            Me(k_2  ,l_2)    = Me(k_2  ,l_2)     +  M_2  *mass_kl(k,l)*I_3*ndp_dx*wp(i);
        end
    end
end
