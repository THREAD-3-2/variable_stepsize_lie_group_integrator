function [g_int_vec] = int_constraints(q, param)
% internal constraints
%
% :param q: current configuration
% :param param: parameters 
%
% :returns: internal constraints

q_tmp = reshape(q,15,param.n_kn_nodes);

g_int = cell(param.n_kn_nodes,1);
%g_int = zeros(6, param.n_kn_nodes);

for i=1:param.n_kn_nodes
    d1 = q_tmp(4:6,i);
    d2 = q_tmp(7:9,i);
    d3 = q_tmp(10:12,i);
    
    g1 = (d1'*d1  -1)/2;
    g2 = (d2'*d2  -1)/2;
    g3 = (d3'*d3  -1)/2;
    g4 = d1'*d2;
    g5 = d3'*d1;
    g6 = d2'*d3;
    
    g_int{i} = [g1;g2;g3;g4;g5;g6];
end

g_int_vec = [g_int{:}];
g_int_vec = reshape(g_int_vec, param.n_kn_nodes*6,1);

end

