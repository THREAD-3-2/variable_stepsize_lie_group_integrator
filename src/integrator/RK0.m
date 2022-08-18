function [out, dout] = RK0(X, q_n, param, fns)
% evaluate residual and tangent for step 1
%
% :param X: configuraiton at step n+1
% :param q_n: configuraiton at step n
% :param param: parameters
% :param fns: symbolic functions
%
% :returns: residual and tangent
q_np1 = X(1:param.n_kn_nodes*15);

R = fns.R_sm(q_n, q_np1);

Kq = fns.Kq_sm(q_n, q_np1);

out = full(R); % R
dout = full(Kq); % K
end