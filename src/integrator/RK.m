function [out, dout] = RK(X, q_nm1, q_n, param, fns)
% evaluate residual and tangent for step n
%
% :param X: configuraiton at step n+1
% :param q_nm1: configuraiton at step n-1
% :param q_n: configuraiton at step n
% :param param: parameters
% :param fns: symbolic functions
%
% :returns: residual and tangent

q_np1 = X(1:param.n_kn_nodes*15);

R = fns.R_vi(q_nm1, q_n, q_np1);

Kq = fns.Kq_vi(q_nm1, q_n, q_np1);

out = full(R); % R
dout = full(Kq); % K

end

