
pminus_fnc=fncs.pminus_fnc;
G_fnc = fncs.G_fnc;

dur = param.N_timesteps;

time = (1:1:dur).*param.timestep;

extEpot = zeros(dur,1);
intEPot = zeros(dur,1);
KE_Hami_red=zeros(dur,1);

for i=1:dur-1
    qn_tmp = res.Q(:,i);
    qnp1_tmp = res.Q(:,i+1);
    
    % momentum p-
    pminus = param.timestep*pminus_fnc(qn_tmp, qnp1_tmp); % vector
    
    % Jacobi G
    Gn = G_fnc(qn_tmp);
    
    % kinetic energy in terms of momentum
    KE_Hami_red(i) = Hami_Energy_red(full(pminus), param, full(Gn));
    
    % external potential energy
    extEpot(i) = param.ext_pot_en_toggle*ext_potEn(qn_tmp,param.q_ref,param);
    
    % internal potential energy
    intEPot(i) =  param.int_pot_en_toggle*int_potEn(qn_tmp,param.q_ref,param);
    
end

total_En = extEpot + intEPot + KE_Hami_red; % total energy

% plot total energy
plot(time(1:end-1),total_En(1:end-1),'r','LineWidth',2)
box on;
% ylim([0.081 0.0825]); 
xlabel('Time / ms')
ylabel('Hamiltonian / mJ')

% plot kinetic enetgy and potential energy
% [AX,H1,H2] = plotyy(time(1:end-1),intEPot(1:end-1),time(1:end-1),KE_Hami_red(1:end-1),'plot');
% legend([H1,H2],{'V_d (\eta=0)','T_d (\eta=0)'})
% set(get(AX(1),'Ylabel'),'String','Potential energy / mJ')
% set(get(AX(2),'Ylabel'),'String','Kinetic energy / mJ')
% xlabel('Time / ms')
% AX(1).YLim = [814e-4 817e-4];
% AX(1).YTick = 0.0814:0.00005:0.0817;
% AX(2).YLim = [-0.0001 0.0004];
% AX(2).YTick = -0.0001:0.0001:0.0004;
% set(H1,'Color','b','LineWidth',2)
% set(H2,'Color','r','LineStyle',':','LineWidth',2)
% set(AX,{'ycolor'},{'b';'r'})
% grid on;
% box on;
% set(gca,'FontSize',15)
% AX(1).FontSize=15;
% AX(2).FontSize=15;
