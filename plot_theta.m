% plot rotation angle of the end node

dur = param.N_timesteps;
theta=zeros(size(dur,2));

for i=1:dur
 j=(param.n_kn_nodes-1)*15+4; % 1-x, 2-y, 3-z, 4-d1x, 5-d1y, 6-d1z, 
 d1=res.Q(j:j+2,i);
 d10=q_ref(j:j+2);
 theta(i)=acosd(dot(d1,d10));
end

time = (1:1:dur).*param.timestep;

hold on;
plot(time, theta,'r','LineWidth',2);
box on;
xlabel('Time / ms')
ylabel('Rotation angle of the end section / o')
hold off;
