function res = exponentialSE3N(sigma)

    N = floor(length(sigma)/6);
    res = zeros(3*N,4);
    for i = 1:N
       res(3*i-2:3*i, :) = expSE3(sigma(6*i-5:6*i));
    end

end