function x = null_space_matrix(q_bm, param)
% null space matrix
%
% :param q_bm: current configuraiton
% :param param: parameters
%
% :returns: null space matrix P

% initialize matrix
dim_q_ng = param.dim_q - 6*param.n_kn_nodes;
x=zeros(param.dim_q,dim_q_ng); % param.dim_q_ng=param.dim_q-num of total constraints

%% ================================================================
% null space matrix, refers to internal constraint: triads is perpendicular

%-----------------------------------------------------------------------------------------------
% first node, total fixed
in=1;
iq_b1_phi__b1_d3(in,:)=[param.iq.bm(in).phi, param.iq.bm(in).d1, param.iq.bm(in).d2, param.iq.bm(in).d3, param.iq.bm(in).v];
iu_b1_wphi__b1_om(in,:)=[param.iu.bm(in).wphi, param.iu.bm(in).om,param.iu.bm(in).ov];

% null space matrix for beam node, Leye.phd thesis p94-eq(5.4.32)
x(iq_b1_phi__b1_d3(in,:),iu_b1_wphi__b1_om(in,:))= zeros(15,9);

%-----------------------------------------------------------------------------------------------
for in =2:param.n_kn_nodes
    iq_b1_phi__b1_d3(in,:)=[param.iq.bm(in).phi, param.iq.bm(in).d1, param.iq.bm(in).d2, param.iq.bm(in).d3, param.iq.bm(in).v];
    iu_b1_wphi__b1_om(in,:)=[param.iu.bm(in).wphi, param.iu.bm(in).om,param.iu.bm(in).ov];
    
    % null space matrix for beam node, Leye.phd thesis p94-eq(5.4.32)
    x(iq_b1_phi__b1_d3(in,:),iu_b1_wphi__b1_om(in,:))=...
        [eye(3),       zeros(3,3),                                         zeros(3,3);
        zeros(3,3),  -hat_vec(q_bm(param.iq.bm(in).d1)),   zeros(3,3);
        zeros(3,3),  -hat_vec(q_bm(param.iq.bm(in).d2)),   zeros(3,3);
        zeros(3,3),  -hat_vec(q_bm(param.iq.bm(in).d3)),   zeros(3,3);
        zeros(3,3),  zeros(3,3),                                         eye(3)];
end % loop over nodes

%-----------------------------------------------------------------------------------------------
% change last eye(3) to zero, if boundary condition imposed on electric dof
for in =2:param.n_kn_nodes
    iq_b1_phi(in,:)=param.iq.bm(in).v;
    iu_b1_wphi(in,:)=param.iu.bm(in).ov;
    switch param.eboun
        case 'end'
            if(in==param.n_kn_nodes) 
                x(iq_b1_phi(in,:),iu_b1_wphi(in,:))=zeros(3,3); % boundary at last node
            end
        otherwise
             % position with elec boundary
            if(rem(in-1,param.nec)==0)
                x(iq_b1_phi(in,:),iu_b1_wphi(in,:))=zeros(3,3);
            end
    end
end
