function [H] = Hami_Energy_red(pminus, param, Gn)
% reduced Kinetic energy: p*invM*p/2
%
H = 0;
p = reshape(pminus,15,param.n_kn_nodes);

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

for i=1:param.n_el_nodes % *global assembly*
    
    % reduced momentum at node i, i+1
    p_i = p(1:9,i);
    p_ip1 = p(1:9,i+1);
    p_e=[p_i; p_ip1]; % 18*1
    
    % reduced constraint Jacobi G,   G_6node*15node
    % G_i =     Gn( 6*(i-1)+1 : 6*i,        15*(i-1)+1 : 15*i);       % 6*15
    % G_i+1 = Gn( 6*i+1      :  6*(i+1), 15*i+1       : 15*(i+1)); % 6*15
    G_iip1 = Gn( 6*(i-1)+1:  6*(i+1), 15*(i-1)+1 : 15*(i+1));  % 12*30
    G_red = [G_iip1(1,1:9)   G_iip1(1,16:24); 
                  G_iip1(2,1:9)   G_iip1(2,16:24);
                  G_iip1(4,1:9)   G_iip1(4,16:24); 
                  G_iip1(7,1:9)   G_iip1(7,16:24); 
                  G_iip1(8,1:9)   G_iip1(8,16:24); 
                  G_iip1(10,1:9) G_iip1(10,16:24)]; % 6*18
    Q = eye(18) - G_red' * inv(G_red*param.Me_bar_inv(:,:,i)*G_red') * G_red * param.Me_bar_inv(:,:,i);
    Qp_e = Q*p_e;
    
    % kinetic energy of element
    T_i = (Qp_e')*param.Me_bar_inv(:,:,i)*(Qp_e)/2;
    
    H = H + T_i; % sum all elemental comtribution *global assembly*
end

end