function plot_FE(q ,param, res, item)
% plot configuration of the beam
%
%

h = param.height; w = param.width;

q = reshape(q, 15, param.n_kn_nodes);
q_0=reshape(param.q_ref, 15, param.n_kn_nodes);


colormap(jet)
hold on
for i_el = 1:param.n_el_nodes
    
    % configuration of node i
    phi_l0 = q_0(1:3,i_el);
    phi_l = q(1:3,i_el);
    d_1_l = q(4:6,i_el);
    d_2_l = q(7:9,i_el);
    v_l     = q(13:15,i_el);
    
    % configuration of node i+1
    phi_r0 = q_0(1:3,i_el+1);
    phi_r = q(1:3,i_el+1);
    d_1_r = q(4:6,i_el+1);
    d_2_r = q(7:9,i_el+1);
    v_r     = q(13:15,i_el+1);
    
    %  vertices of corss section i
    c1 =phi_l + h(i_el)/2*d_2_l + w(i_el)/2*d_1_l;
    c2 =phi_l - h(i_el)/2*d_2_l + w(i_el)/2*d_1_l;
    c3 =phi_l - h(i_el)/2*d_2_l - w(i_el)/2*d_1_l;
    c4 =phi_l + h(i_el)/2*d_2_l - w(i_el)/2*d_1_l;
    
    %  vertices of corss section i+1
    c5 =phi_r + h(i_el+1)/2*d_2_r + w(i_el+1)/2*d_1_r;
    c6 =phi_r - h(i_el+1)/2*d_2_r + w(i_el+1)/2*d_1_r;
    c7 =phi_r - h(i_el+1)/2*d_2_r - w(i_el+1)/2*d_1_r;
    c8 =phi_r + h(i_el+1)/2*d_2_r - w(i_el+1)/2*d_1_r;
    
    v = [c1 c2 c3 c4 c5 c6 c7 c8]';
    
    % face
     f = [1 4 3 2;...
          5 6 7 8;...
          1 2 6 5;...
          3 4 8 7;...
          2 3 7 6;...
          1 5 8 4];
    % plot
    switch item
        case 'u'
            % displacement z at nodes
            u_l = phi_l-phi_l0;
            u_r = phi_r-phi_r0;
            
            %  vertices of corss section i
            col_1 = u_l(3) + h(i_el)/2*d_2_l(3) + w(i_el)/2*d_1_l(3);
            col_2 = u_l(3) - h(i_el)/2*d_2_l(3) +  w(i_el)/2*d_1_l(3);
            col_3 = u_l(3) -  h(i_el)/2*d_2_l(3) -  w(i_el)/2*d_1_l(3);
            col_4 = u_l(3) +  h(i_el)/2*d_2_l(3) - w(i_el)/2*d_1_l(3);
            
            %  vertices of corss section i+1
            col_5 = u_r(3) + h(i_el+1)/2*d_2_r(3) + w(i_el+1)/2*d_1_r(3);
            col_6 = u_r(3) - h(i_el+1)/2*d_2_r(3) +  w(i_el+1)/2*d_1_r(3);
            col_7 = u_r(3) -  h(i_el+1)/2*d_2_r(3) -  w(i_el+1)/2*d_1_r(3);
            col_8 = u_r(3) +  h(i_el+1)/2*d_2_r(3) - w(i_el+1)/2*d_1_r(3);
            
            col = [col_1;col_2;col_3;col_4;col_5;col_6;col_7;col_8]; % color on vertices
            
            % FaceColor: flat(color per surface) or interp(color per vertex); , 'EdgeColor', 'red',
            patch('Faces',f, 'Vertices',v,'FaceVertexCData',col,  'FaceColor','interp');
            
        case 'v'
            % electrical potential on the vertices
            %  vertices of corss section i
            col_1 = v_l(1) + h(i_el)/2*v_l(3) + w(i_el)/2*v_l(2);
            col_2 = v_l(1) - h(i_el)/2*v_l(3) +  w(i_el)/2*v_l(2);
            col_3 = v_l(1) -  h(i_el)/2*v_l(3) -  w(i_el)/2*v_l(2);
            col_4 = v_l(1) +  h(i_el)/2*v_l(3) - w(i_el)/2*v_l(2);
            
            %  vertices of corss section i+1
            col_5 = v_r(1) + h(i_el+1)/2*v_r(3) + w(i_el+1)/2*v_r(2);
            col_6 = v_r(1) -  h(i_el+1)/2*v_r(3) + w(i_el+1)/2*v_r(2);
            col_7 = v_r(1) -  h(i_el+1)/2*v_r(3) -  w(i_el+1)/2*v_r(2);
            col_8 = v_r(1) + h(i_el+1)/2*v_r(3) -  w(i_el+1)/2*v_r(2);
            
            col = [col_1;col_2;col_3;col_4;col_5;col_6;col_7;col_8]; % color on vertices
            
            % FaceColor: flat(color per surface) or interp(color per vertex); , 'EdgeColor', 'red',
            patch('Faces',f, 'Vertices',v,'FaceVertexCData',col,  'FaceColor','interp'); 
        case 'm'
            % bending and torsion stress
            col_1 = res.stress_m(3*i_el,end); % 3*i_el: m-z;  3*i_el-1: m-y; 3*i_el-2: m -x
            
            col = [col_1;col_1;col_1;col_1;col_1;col_1]; % color on surfaces
            
            patch('Faces',f, 'Vertices',v,'FaceVertexCData',col,  'FaceColor','flat');
        case 'n'
            % shear and elongation
            col_1 = res.stress_n(3*i_el-2,end);
            
            col = [col_1;col_1;col_1;col_1;col_1;col_1]; % color on surfaces
            
            patch('Faces',f, 'Vertices',v,'FaceVertexCData',col,  'FaceColor','flat'); 
        otherwise
            error('no approximation selected')
    end
end

hold off

% set axis
phi_i_x(1,:) = q(1,:);
phi_i_y(1,:) = q(2,:);
phi_i_z(1,:) = q(3,:);

xlim_plus = 1.1*max(phi_i_x,[],'all')+max(param.height);
ylim_plus = 1.3*max(phi_i_y,[],'all')+max(param.width);
zlim_plus = 1.3*max(phi_i_z,[],'all');
xlim_minus = 1.1*min(phi_i_x,[],'all')-max(param.height);
ylim_minus = 1.3*min(phi_i_y,[],'all')-max(param.width);
zlim_minus = 1.3*min(phi_i_z,[],'all');
lim = [xlim_minus; xlim_plus; ylim_minus; ylim_plus; zlim_minus; zlim_plus];
axis equal
axis(lim)
xlabel('x')
ylabel('y')
zlabel('z')
grid on
% campos([-190, 150, 140])
view(45,23)

% set color bar
colorbar
clb = colorbar;
switch item
    case 'u'
        ylabel(clb, 'Displacement_z / mm');
    case 'v'
        ylabel(clb, 'Electric potential / V');
    case 'm'
        ylabel(clb, 'Bending and torsion');
    case 'n'
        ylabel(clb, 'Shear force n_x');
end
caxis('auto');

end
