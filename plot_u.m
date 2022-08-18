% plot position/displacement of the end node (0,0,0.1)

dur = param.N_timesteps;
q_end=zeros(size(dur,2));

for i=1:dur
 j=(param.n_kn_nodes-1)*15+3; % 1-x,2-y,3-z,
 q_end(i)=res.Q(j,i);
end

hold on;
plot(param.time, q_end,'r','LineWidth',2);
box on;
xlabel('Time / ms')
ylabel('Position of the end node / mm')
hold off;

xlim([0 0.1]);
ylim([0.09 0.101]);