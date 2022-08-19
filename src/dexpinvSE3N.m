function res = dexpinvSE3N(sigma, input)
% dexpinv map on SE(3)^N
%
% :param sigma: 
% :param input:

% :returns: 

    N = floor(length(sigma)/6);
    res = zeros(6*N,1);
    for i = 1:N
       res(6*i-5:6*i)= dexpinvSE3(sigma(6*i-5:6*i),input(6*i-5:6*i));
    end

end
