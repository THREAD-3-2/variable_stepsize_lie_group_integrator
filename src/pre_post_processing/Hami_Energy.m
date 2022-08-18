function [H] = Hami_Energy(pminus, param)
% Kinetic energy: p*invM*p/2
%
%
H = 0;
p = reshape(pminus,15,param.n_kn_nodes);

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

for i=1:param.n_el_nodes % *global assembly*
    
    % reduced momentum at node i, i+1
    p_i = p(1:9,i);
    p_ip1 = p(1:9,i+1);
    p_e=[p_i; p_ip1];
    
    % kinetic energy of element
    T_i = (p_e')*param.Me_bar_inv(:,:,i)*(p_e)/2;
    
    H = H + T_i; % sum all elemental comtribution *global assembly*
end

end