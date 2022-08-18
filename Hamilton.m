% plot discrete Hamiltonian

dur = param.N_timesteps;
extEpot = zeros(dur,1);
intEPot = zeros(dur,1);
KE_Hami=zeros(dur,1);

for i=1:dur-1
    qn_tmp = res.Q(:,i);
    qnp1_tmp = res.Q(:,i+1);
    
    % momentum p-
    pminus = param.timestep*fns.pminus(qn_tmp, qnp1_tmp); % vector
    
    % kinetic energy in terms of momentum
    KE_Hami(i) = Hami_Energy(full(pminus), param);
    
    % external potential energy
    extEpot(i) = 0*ext_potEn(qn_tmp,param.q_ref,param);
    
    % internal potential energy
    intEPot(i) =  int_potEn(qn_tmp,param.q_ref,param);
end

total_En = extEpot + intEPot + KE_Hami;% total energy

% plot total energy
plot(param.time(1:end-1),total_En(1:end-1),'r','LineWidth',2)
box on;
xlabel('Time / ms')
ylabel('H_d / mJ')
xlim([0 0.1]);
ylim([0.081 0.0825]);
