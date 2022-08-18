function [W_int] = int_potEn(q, q_ref, param)
% internal potential energy
%
% :param q: current configuraiton
% :param q_ref: reference configuraiton
% :param param: parameters
%
% :returns: internal potential energy

W_int = 0;

q_ref_tmp = reshape(q_ref, 15, param.n_kn_nodes);
q_tmp = reshape(q, 15, param.n_kn_nodes);

for i=1:param.n_el_nodes % *global assembly*
    
    e_ref_i = [q_ref_tmp(:,i) q_ref_tmp(:,i+1)]; % two nodes q1_0, q2_0
    
    e_i = [q_tmp(:,i) q_tmp(:,i+1)]; % two nodes: q1_n, q_2n
    
    W_int_i = int_potEn_bm(e_ref_i, e_i, param); % potential energy of element i
    
    W_int = W_int + W_int_i; % sum all elemental comtribution *global assembly*
    
end

end

