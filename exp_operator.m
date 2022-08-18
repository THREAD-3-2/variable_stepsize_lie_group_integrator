function Lambda=exp_operator(Theta)
% exponential map
%
theta=norm(Theta);
if theta == 0
    Lambda=eye(3);
else
    c_1=sin(theta); c_2=cos(theta); n=Theta/theta;
    Lambda=c_2*eye(3)+(1-c_2)*n*n'+c_1*hat_vec(n);
end