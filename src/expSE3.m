function A = expSE3(input)

    %we take as an input the element of SE(3) representable by [hat(u),v]
    %so u and v are column 3-vectors
    
    u = input(1:3);
    v = input(4:6);
    theta = norm(u,2);
    
    tol = 1e-16;  
    
    if theta>tol
        A = sin(theta)/theta;
        B = (1-cos(theta))/(theta^2);
        C = (1-A)/(theta^2);
        V = eye(3) + B*hat(u) + C * hat(u) * hat(u);
        A = [expRodrigues(u), V*v];
    elseif theta==0
        A = [expRodrigues(u), v];
    else
        Blow = 0.5-theta^2/24 + theta^4/720 - theta^6/40320;
        Clow = (1/6-theta^2/120+theta^4/5040-theta^6/362880);
        Vlow = eye(3) + Blow*hat(u) + Clow * hat(u) * hat(u);
        A = [expRodrigues(u), Vlow*v];
    
    end
    
end