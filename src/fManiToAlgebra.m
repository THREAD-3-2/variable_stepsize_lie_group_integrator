function F = fManiToAlgebra(q,w,L,m)

%     R = assembleR(q,L,m); %Creates the matrix R multiplying w'
%     
%     Func = assembleF(q,w,m,L); %Creates the right hand site Rw' = Func
    
    vec = @(v,i) getVec(v,i);
    
    N = length(m); %Number clcof connected pendulums
    
%     V = R\Func; %Finds the right hand side of the equation w' = V
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