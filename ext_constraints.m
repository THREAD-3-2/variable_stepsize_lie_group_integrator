function [g_ext] = ext_constraints(q, q_ref, param)
% external constraints
%
% :param q: current configuraiton
% :param q_ref: reference configuraiton
% :param param: parameters
%
% :returns: external constraints

q_ref_tmp = reshape(q_ref,15,param.n_kn_nodes);
q_tmp = reshape(q,15,param.n_kn_nodes);

g_ext = [ q_tmp(1:3,1) - q_ref_tmp(1:3,1); %fixing one node
    (q_tmp(4:6,1)'*q_ref_tmp(7:9,1));
    (q_tmp(7:9,1)'*q_ref_tmp(10:12,1));
    (q_tmp(10:12,1)'*q_ref_tmp(4:6,1)); %preventing node rotation
    q_tmp(13,1)-0; % 0 kV
    q_tmp(14,1)-0; % 0 kV/mm
    q_tmp(15,1)-0; % 0 kV/mm
    %                        q_bm1(13,2)-param.ele_pot;
    q_tmp(13,end)-param.ele_pot(1); % 0 kV at beam node
    q_tmp(14,end)-param.ele_pot(2); % 0 kV/mm
    q_tmp(15,end)-param.ele_pot(3)]; % 0 kV/mm
end