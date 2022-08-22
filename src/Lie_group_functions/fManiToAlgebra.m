function F = fManiToAlgebra(q,w,L,m)
% RHS of the system
%
% :param q: 
% :param w: 
% :param L: 
% :param m: 
%
% :returns: vector field of the system

    vec = @(v,i) getVec(v,i);
    
    N = length(m); % Number of connected pendulums
    
    z = [q;w];
    V = FuncW(z,L,m);
    A = @(i) hat(vec(q,i))*vec(V,i);
    F = zeros(6*N,1);
   
    for i = 1:2*N
        if mod(i,2)==0
            F(3*i-2:3*i) =  A(i/2);
        else
            F(3*i-2:3*i) = vec(w,(i+1)/2) ;
        end
    end
    
end
