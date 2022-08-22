function w = extractw(v)
    
    %from the input vector v = [q1,w1,q2,w2,...,qN,wN] it gives as an
    %output just the w components : w = [w1,w2,...,wN]

    N = length(v)/6; %Number of connected pendulums
    w = zeros(3*N,1);
    for i = 1:N
        w(3*i-2:3*i) = v(6*i-2:6*i);
    end
end