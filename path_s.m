function [q] = path_s(L, param)
% generates initial configutaiton q for straight beam with length L
%

%% beam geometry (length, width, height in [mm])

syms p_x(s)
syms p_y(s)
syms p_z(s)

%% path definition
e_x = [1; 0; 0]; e_y = [0; 1; 0]; e_z = [0; 0; 1];
p_x(s) = 0;
p_y(s) = 0;
p_z(s) = s;

%differentiation for directors
 dpds_x = diff(p_x,s,1); %tangent -> d_3
ddpds_x = diff(p_x,s,2); %normal -> d_2
 dpds_y = diff(p_y,s,1); %tangent -> d_3
ddpds_y = diff(p_y,s,2); %normal -> d_2
 dpds_z = diff(p_z,s,1); %tangent -> d_3
ddpds_z = diff(p_z,s,2); %normal -> d_2

%% path discretization
discretization = linspace(0,1,param.n_kn_nodes);
phi_disc = double([p_x(discretization);p_y(discretization);p_z(discretization)])*L;

dpds_disc = double([ dpds_x(discretization); dpds_y(discretization); dpds_z(discretization)]);
ddpds_disc = double([ddpds_x(discretization);ddpds_y(discretization);ddpds_z(discretization)]);

phi_e = [discretization*0; discretization*0; discretization*0];% initial electrica potential

%% directors
% director definition, d_3 = tangent, d_2 || to e_z, d_1 || e_y
d_3 = dpds_disc;
d_1 = -cross(d_3,e_y.*ones(1,param.n_kn_nodes));
d_2 = cross(d_3,d_1);

% scaling
d_3 = d_3./sqrt( d_3(1,:).^2 +  d_3(2,:).^2 +  d_3(3,:).^2);
d_2 = d_2./sqrt( d_2(1,:).^2 +  d_2(2,:).^2 +  d_2(3,:).^2);
d_1 = d_1./sqrt( d_1(1,:).^2 +  d_1(2,:).^2 +  d_1(3,:).^2);

q = [phi_disc; d_1; d_2; d_3; phi_e]; % *last component is electrical potential*
% [phi_1 phi_2 ...  phi_n]
% [ d1_1  d1_2 ...  d1_n]
% [ d2_1  d2_2 ...  d2_n]
% [ d3_1  d3_2 ...  d3_n]
% [ eo_1  eo_2 ...  eo_n]
% [ ex_1  ex_2 ...  ex_n]
% [ ey_1  ey_2 ...  ey_n]

%% reshaping
q = reshape(q, 1, []);
% disp(q)
end

