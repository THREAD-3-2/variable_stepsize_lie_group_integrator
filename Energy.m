% plot int pot energy and kinetic energy

dur = param.N_timesteps;

extEpot = zeros(dur,1);
intEPot = zeros(dur,1);
KinEn=zeros(dur,1);

for i=1:dur
    q_tmp = res.Q(:,i);
    
    % internal potential energy
    intEPot(i) =  int_potEn(q_tmp,param.q_ref,param);
    
    if i>1
        q_dot = (q_tmp - res.Q(:,i-1))/param.timestep; % finite difference approaximation of q_dot
        
        KinEn(i) = kinEnergy_new(q_dot, param);
    end
    
end
total_En = intEPot + KinEn;% total energy

time = (1:1:dur).*param.timestep;
box on;

%plot total energy
plot(time(1:end-1),total_En(1:end-1),'b','LineWidth',2)
% ylim([0.0815 0.0816]); % for total energy
xlabel('Time / ms')
ylabel('Energy / mJ')
