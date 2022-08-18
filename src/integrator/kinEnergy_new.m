function [E_kin] = kinEnergy_new(q_dot, param)
% kinetic energy
%
% :param q_dot: velocity
% :param param: parameters
%
% :returns: kinetic energy

E_kin = 0;

q_dot_tmp = reshape(q_dot,15,param.n_kn_nodes);

for i=1:param.n_el_nodes % *global assembly*
    
    % velocity at node i, i+1
    q_i_dot_n = q_dot_tmp(:,i);
    q_ip1_dot_n = q_dot_tmp(:,i+1);
    q_dot = [q_i_dot_n; q_ip1_dot_n];
    
    % kinetic energy of element
    E_kin_i = (q_dot')*param.Me(:,:,i)*(q_dot)/2;
    
    E_kin = E_kin + E_kin_i; % sum all elemental comtribution *global assembly*
end

end

