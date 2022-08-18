function plot_FE_mov(Q , param)
% creates a movie from Q and colors the beam if strain results are
% provided, saved to code folder under sim_vid

% INPUTS:
% Q: beam config over time
% param: struct with parameters
% node_col: avaraged coloring value per node for each timestep

% chose movie FPS
fps = 30;

h = param.height; w = param.width;

dur = param.N_timesteps;

%% axis scaling
fig = figure('Color','white','units','normalized','outerposition',[0 0 1 1]);
vid = VideoWriter('Sim_vid', 'MPEG-4');
open(vid);

% t_i = 1:floor(dur/(param.totaltime*fps)):dur;
for t_i = 1:dur
    % clearing frame
    if t_i<dur
        clf('reset')
    end
    
    colormap(jet)
    q = Q(:,t_i);
    q = reshape(q, 15, param.n_kn_nodes);
    q_0=reshape(param.q_ref, 15, param.n_kn_nodes);
    hold on
    % set axis
    phi_i_x(1,:) = q_0(1,:);
    phi_i_y(1,:) = q_0(2,:);
    phi_i_z(1,:) = q_0(3,:);
    
    xlim_plus = 1.1*max(phi_i_x,[],'all')+max(param.height);
    ylim_plus = 1.3*max(phi_i_y,[],'all')+max(param.width);
    zlim_plus = 1.3*max(phi_i_z,[],'all');
    xlim_minus = 1.1*min(phi_i_x,[],'all')-max(param.height);
    ylim_minus = 1.3*min(phi_i_y,[],'all')-max(param.width);
    zlim_minus = 1.3*min(phi_i_z,[],'all');
    lim = [xlim_minus; xlim_plus; ylim_minus; ylim_plus; zlim_minus; zlim_plus];
    axis(lim)
    axis equal;
    xlabel('x')
    ylabel('y')
    zlabel('z')
    grid on
    %     campos([150,  -190,  140])
    view(45,23)
    
    % set color bar
    colorbar
    caxis('auto');
    
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
        c2 =phi_l + h(i_el)/2*d_2_l - w(i_el)/2*d_1_l;
        c3 =phi_l - h(i_el)/2*d_2_l - w(i_el)/2*d_1_l;
        c4 =phi_l - h(i_el)/2*d_2_l + w(i_el)/2*d_1_l;
        
        %  vertices of corss section i+1
        c5 =phi_r + h(i_el+1)/2*d_2_r + w(i_el+1)/2*d_1_r;
        c6 =phi_r + h(i_el+1)/2*d_2_r - w(i_el+1)/2*d_1_r;
        c7 =phi_r - h(i_el+1)/2*d_2_r - w(i_el+1)/2*d_1_r;
        c8 =phi_r - h(i_el+1)/2*d_2_r + w(i_el+1)/2*d_1_r;
        
        v = [c1 c2 c3 c4 c5 c6 c7 c8]';
        
        % face
        f = [1 2 3 4;...
            5 6 7 8;...
            1 2 6 5;...
            3 4 8 7;...
            2 3 7 6;...
            1 4 8 5];
        
        % electrical potential on the vertices
        %  vertices of corss section i
        col_1 = v_l(1) + h(i_el)/2*v_l(3) + w(i_el)/2*v_l(2);
        col_2 = v_l(1) - h(i_el)/2*v_l(3) +  w(i_el)/2*v_l(2);
        col_3 = v_l(1) -  h(i_el)/2*v_l(3) -  w(i_el)/2*v_l(2);
        col_4 = v_l(1) +  h(i_el)/2*v_l(3) - w(i_el)/2*v_l(2);
        
        %  vertices of corss section i+1
        col_5 = v_r(1) + h(i_el+1)/2*v_r(3) + w(i_el+1)/2*v_r(2);
        col_6 = v_r(1) - h(i_el+1)/2*v_r(3) +  w(i_el+1)/2*v_r(2);
        col_7 = v_r(1) -  h(i_el+1)/2*v_r(3) -  w(i_el+1)/2*v_r(2);
        col_8 = v_r(1) +  h(i_el+1)/2*v_r(3) - w(i_el+1)/2*v_r(2);
        
        col = [col_1;col_2;col_3;col_4;col_5;col_6;col_7;col_8]; % color on vertices
        % ++++++++++++++++++++++++++++++++++++++++++++++++++++
        patch('Faces',f, 'Vertices',v,'FaceVertexCData',col,  'FaceColor','interp'); %alpha(0.1)
    end
    pause(param.timestep)
    %     pause(1)
    
    % writing frame to vid
    frame = getframe(fig);
    writeVideo(vid,frame);
    
end
hold off
close(vid)

end

