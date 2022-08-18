function [g] = constraints(q, q_ref, param)
% constraints
%
% :param q: current configuraiton
% :param q_ref: reference configuraiton
% :param param: parameters
%
% :returns: constraints

% g = [int_constraints(q, param); ext_constraints(q, q_ref, param)];
g = [int_constraints(q, param)];

end