function [err, TT, Y] = variableStepComparison(vecField, action, z0, T, tol) 

    chunk = 500;
    Z = zeros(length(z0), chunk);
    TT = zeros(1, chunk);
    Y = Z;
    Y(:, 1) = z0;

    a = 1/4;
    theta = 0.85;
    i = 1;
    h = T/499;
    rejected = 0;

    while TT(i) < T - 5 * eps
        err = tol + 1;
        while err > tol    
            [z, err] = variableRKMK45(vecField, action, Y(:, i), h);
            accepted = (err < tol);
            if accepted
                i = i + 1;
                Y(:, i) = z;
                TT(i) = TT(i - 1) + h;
            else
                rejected = rejected + 1;
            end
            h = min(theta * (tol/err)^a * h, T - TT(i));
 
        end
        if mod(i, chunk) == 0
            Y = [Y, Z];
            TT = [TT zeros(1, chunk)];
        end
    end
    Y = Y(:, 1 : i);
    TT = TT(: , 1 : i); 
    
end