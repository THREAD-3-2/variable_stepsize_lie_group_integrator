function [W_ext] = ext_potEn(q, q_ref, param)
% external potential energy
%
% :param q: current configuraiton
% :param _ref: reference configuraiton
%
% :returns: external potential energy

% gravity vector,
g_vec = param.g_dir*param.g; %g in negative z-axis

W_ext = 0;

q_tmp = reshape(q, 15, param.n_kn_nodes);

for i=1:param.n_el_nodes % *global assembly*
    
    e_i = [q_tmp(:,i) q_tmp(:,i+1)]; % two nodes: q1_n, q_2n
    
    h1=[0 0 e_i(3,1) 0 0 0 0 0 0 0 0 0 0 0 0]';
    h2=[0 0 e_i(3,2) 0 0 0 0 0 0 0 0 0 0 0 0]';
    h = [h1; h2];
    
    W_ext_i = [g_vec g_vec]*param.Me(:,:,i)*h; % potential energy of element i
    
    W_ext = W_ext + W_ext_i; % sum all elemental comtribution *global assembly*
    
end

end