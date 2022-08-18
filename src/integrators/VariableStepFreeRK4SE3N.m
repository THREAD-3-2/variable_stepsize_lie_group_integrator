function [y,err] = VariableStepFreeRK4SE3N(f,action,h,p)
        
    gAc = @(g,x) action(exponentialSE3N(g),x); 
    
    Y1 = p;
    k1 = h*f(Y1);
    
    Y2 = gAc(k1/2,p);
    k2 = h*f(Y2);
    
    Y3 = gAc(k2/2,p);
    k3 = h*f(Y3);
    
    Y3bar = gAc(3/4 * k2,p);
    k3bar = h*f(Y3bar);
    yHat = gAc(1/9*(-k1+3*k2+4*k3bar),gAc(1/3*k1,p));
    
    Y4 = gAc(k3-k1/2,Y2);
    k4 = h*f(Y4);
    
    yHalf = gAc(1/12*(3*k1+2*k2+2*k3-k4),p);
    y = gAc(1/12*(-k1+2*k2+2*k3+3*k4),yHalf); 
    
    err = norm(y- yHat,2);

end