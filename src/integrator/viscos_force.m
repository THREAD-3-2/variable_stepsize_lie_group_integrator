function [fext] = viscos_force(q_n, q_np1, q_ref, param)
% viscosity force vector as external forces
%
% :param q_n: configuraiton at step n
% :param q_np1: configuraiton at step n+1
% :param q_ref: reference configuraiton
% :param param: parameters
%
% :returns: viscosity force vector

% temperal descrtize velocity
q_dot = (q_np1 - q_n)/param.timestep;

% midpoint rule for q
q = (q_n + q_np1)/2;

q_ref_tmp = reshape(q_ref,15,param.n_kn_nodes); % 13*nodes
q_tmp = reshape(q,15,param.n_kn_nodes);
q_dot_tmp = reshape(q_dot,15,param.n_kn_nodes);

% initialize the nodal viscos force vector
I_dFdqnP=zeros(size(q_n)); % 78*1
fext = zeros(size(q_n)); % 78*1

%%
for j=1:param.n_el_nodes % loop all elements/elemental assembly
    
    % nodes of element j
    e_ref_j = [q_ref_tmp(:,j) q_ref_tmp(:,j+1)]; % two nodes q1_0, q2_0 related to element j, 13*2
    e_j = [q_tmp(:,j) q_tmp(:,j+1)]; % two nodes: q1_n, q_2n related to element j
    e_dot_j = [q_dot_tmp(:,j) q_dot_tmp(:,j+1)];
    
    % integration over cross section
    X_node=[-param.b/2 -param.b/2; param.b/2 -param.b/2; param.b/2 param.b/2; -param.b/2 param.b/2]; % nodes of the cross section
    alpha=[-1/sqrt(3) 1/sqrt(3) -1/sqrt(3) 1/sqrt(3)];
    beta =[-1/sqrt(3) -1/sqrt(3) 1/sqrt(3) 1/sqrt(3)];
    wht=1;
    
    % integration over cross section
    for pi=1:4
        Na(:,1)=1/4*(1-alpha(pi))*(1-beta(pi));
        Na(:,2)=1/4*(1+alpha(pi))*(1-beta(pi));
        Na(:,3)=1/4*(1+alpha(pi))*(1+beta(pi));
        Na(:,4)=1/4*(1-alpha(pi))*(1+beta(pi));
        XY=Na*X_node; % paeameterized X1, X2 on cross section
        dNa=[-1/4*(1- beta(pi)) 1/4*(1-beta(pi)) 1/4*(1+beta(pi)) -1/4*(1+beta(pi));
            -1/4*(1-alpha(pi)) -1/4*(1+alpha(pi)) 1/4*(1+alpha(pi)) 1/4*(1-alpha(pi))];
        Ja1=det((dNa*X_node)'); % Jacobian dX/dxi
        
        % viscos stress P and deformation gradient F in element j - function
        [P_e_vis, F_e, Ja2]= viscos_force_PF(e_ref_j, e_j, e_dot_j, XY, param); % Ja2=ds/dxi
        
        % rewrite 3*3 tensor to 9*1 vector
        F = reshape(F_e,[9,1]);
        P_vis = reshape(P_e_vis,[9,1]);
        
        % dF/dqn
        dF_dqn = jacobian(F, q_n); % 9*78
        
        % double contraction dF/dqn : P
        dFdqnP=(dF_dqn(1,:).*0)';
        for jj=1:size(q_n,1)
            for ii=1:9
                dFdqnP(jj) = dFdqnP(jj) + dF_dqn(ii,jj)*P_vis(ii);  % 1*78 vector
            end
        end
        
        % Gauss integration over cross section
        I_dFdqnP = I_dFdqnP + dFdqnP*Ja1*wht;
    end
    
    % Gauss integration over axis s, only one integration point per element
    wp=2;
    fext = fext+ I_dFdqnP*Ja2*wp; % 1*78,    assemble all element by adding up
end

end
