function [L] = LagrangeEqu(q_n, q_np1, q_ref, param)
% discrete Lagrange
%
% :param q_n: configuraiton at step n
% :param q_np1: configuraiton at step n+1
% :param q_ref: configuraiton in reference configuraiton
% :param param: parameters
%
% :returns: discrete Lagrange

% finite diff for q_dot, eq.(72)
q_dot = (q_np1 - q_n)./param.timestep;

% midpoint rule for q, eq.(72)
q = (q_n + q_np1)/2;

%L = T - V_int - V_ext, eq.(66)
L = kinEnergy_new(q_dot, param) - int_potEn(q, q_ref, param) - 0*ext_potEn(q, q_ref, param);
end

