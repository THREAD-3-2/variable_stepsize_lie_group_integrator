function vec = FuncW(z,L,m)
% -------
%
% :param z: 
% :param L: 
% :param m: 
%
% :returns: 

    %This function is used to integrate with ODE45, so the input z is of
    %the form z = [q1,q2,...,qP,w1,w2,...,wP]
    %Builds the part of the vector field for the \dot{w}_i, so
    %R(q)^{-1}*right hand side of the ODE, defined by assembleF.

    q = z(1:length(z)/2);
    w = z(length(z)/2+1:end);
    
    R = assembleR(q,L,m);
    F = assembleF(q,w,m,L);
    
    vec = R\F;
end
