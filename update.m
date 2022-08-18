function x=update(dx, q, param)
% nodal reparametrization
%
% :param dx: generalized configuration
% :param q: redandent configuraiton
% :param param: parameters
%
% :returns: updated redandent configuraiton

x = zeros(15, param.n_kn_nodes);

q = reshape(q, 15, param.n_kn_nodes);

dx = reshape(dx, 9, param.n_kn_nodes);

% !!! here both the internal and external constraints are fullfilled !!!

%% ===============================================================
% first beam node: phi, di and v are fixed in the initial position

x(1:3,1) = q(1:3,1);
x(4:6,1) = q(4:6,1);
x(7:9,1) = q(7:9,1);
x(10:12,1) = q(10:12,1);
switch param.eboun
    case 'end'
        phi0 = 2e4;
        x(13:15,1) = [0, 0, 0]'; % the electric boundary condition: zero
    case 'torsion'
        theta = 0; phi0 = 2e3; alpha = 5e4; beta = 5e4;
        x(13:15,1) = [phi0, alpha, beta]';
    case 'shear'
        phi0 = 2e2; alpha = 2e4;
        x(13:15,1) = [phi0, alpha, 0]';
    case 'bending'
        x(13:15,1) = [0, 0, 0]';
end
%% ===============================================================
% Leye 08a partIII eq(37),
% nodal reparametrisation for mechanical dof: q_n+1 = F(q_n,u), corresponding to P_int
for in = 2:param.n_kn_nodes
    
    dphi = dx(1:3, in); % du
    dtheta = dx(4:6, in);% d_theta/omega
    
    dRot = exp_operator(dtheta);
    
    x(1:3,in) = q(1:3,in) + dphi;
    x(4:6,in) = dRot*q(4:6,in);
    x(7:9,in) = dRot*q(7:9,in);
    x(10:12,in) = dRot*q(10:12,in);
end

%% ===============================================================
% nodal reparametrisation for electric dof:
for in = 2:param.n_kn_nodes
    dv = dx(7:9,in); % d_v
    switch param.eboun
        case 'end'
            % uniaxial contraction, 1st node=[0,0,0], 6 nodes, L0.1,b0.02
            if (in<param.n_kn_nodes)
                x(13:15,in) = q(13:15,in) + dv; % no boun on internal nodes
            else
                x(13:15,in) = [phi0, 0, 0]'; % contraction [2e4, 0, 0]
            end
        case 'torsion'
            % pure torsion, 1st node=[2e3, 3e4, 3e4], 20, 40, 80 nodes,L0.05
            if(rem(in-1,param.nec) == 0)
                theta = theta + 10*pi/param.ncell;
                x(13,in) = phi0;
                x(14,in) = alpha*(cos(theta) - sin(theta));
                x(15,in) = beta*(cos(theta) + sin(theta));
            else
                x(13:15,in) = q(13:15,in) + dv; % no boun on internal nodes
            end
        case 'shear'
            % pure shear, 1st node = [2e2, 2e4, 0], 40 nodes, L0.05b0.005
            if(rem(in-1,param.nec) == 0)
                phi0 = phi0+8e3/param.ncell;
                x(13:15,in) = [phi0, alpha, 0]';
            else
                x(13:15,in) = q(13:15,in) + dv; % no boun on internal nodes
            end
        case 'bending'
            % pure beanding, 1st node=[0, 0, 0]', 40 nodes
            if(rem(in-1,param.nec) == 0)
                if rem((in-1)/param.nec,2) == 0
                    x(13:15,in) = [0, 0, 0]';
                else
                    x(13:15,in) = [2e2, 2e4, 0]';
                end
            else
                x(13:15,in) = q(13:15,in) + dv; % no boun on internal nodes
            end
    end
end

x = reshape(x, [], 1);
end
