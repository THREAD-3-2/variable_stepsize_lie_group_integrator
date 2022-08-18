%This code contains the numerical experiments with the step adaptive RKMK method based on Dormand-Prince pair (5,4)

clc;
clear all;
close all;

%% Setting the parameters 

% prompt = 'How many 3D pendulums do you want to connect?\n\n';
% P = input(prompt);
Prange = 2:2:20;
%Prange = 15

accVar = zeros(length(Prange),1);
accConst = accVar;
index = 1;
tol = 1e-6;

Lref = 5;

steps = zeros(length(Prange),1);


for P = Prange
    
    L = rand(P,1)+0.5; %Set random lengths
    m = rand(P,1)+0.5; %Set random masses
    %L = 0*L + 1;
    L = 0*L + Lref/P; 
    m = 0*m + 1;
    
    P
    
%     [q0,w0,z0] = initializeSE3N(P);
    %[q0,w0,z0] = initializeSE3N_largeVariation(P);
    [q0,w0,z0] = initializeStat(P);
%     [q0,w0,z0] = initializeSE3N_smallVariation(P);
    
    Energy = @(q,w) 0.5*w'*assembleR(q,L,m)*w + potential(q,L,m);

    disp("Energy of this initial condition: "+num2str(Energy(q0,w0)));

    t0 = 0;
    T = 3; 
    N = 1000; 
    time = linspace(t0,T,N); 
    dt = time(2)-time(1);

    getq = @(v) extractq(v);
    getw = @(v) extractw(v);

    f = @(v) fManiToAlgebra(getq(v),getw(v),L,m); 
    action = @(B,input) actionSE3N(B,input); 
    vecField = @(sigma,p) dexpinvSE3N(sigma,f(action(exponentialSE3N(sigma),p)));

    z = z0;
    qC = zeros(3*P,N);
    pC = qC;

    Len = zeros(3*P,1);
    for i = 1:P
        Len(3*i-2:3*i) = L(i)*ones(3,1);
    end
    Mat = diag(Len);
    if P>1
        for i = 3:3:3*(P-1)
            Mat = Mat + diag(Len(1:3*P-i),-i);
        end
    end
    qC(:,1) = q0;
    pC(:,1) = Mat*q0;


    %% COMPARISON BETWEEN CF34 AND CF4, i.e. VARIABLE STEPSIZE AGAINST CONSTANT ONE

    chunk=400;
    Z=zeros(length(z0),chunk);
    TT=zeros(1,chunk);
    Y=Z;
    Y(:,1)=z0;

    a = 1/4;
    theta = 0.85;
    i = 1;
    h = T/100;
    rejected = 0;
    ctr = 0;
    % while TT(i)<T-5*eps
    %     err = tol + 1;
    %     while err>tol    
    %         [z,err] = VariableStepFreeRK4SE3N(f,action,h,Y(:,i));
    %         accepted = (err<tol);
    %         if accepted
    %             i = i + 1;
    %             Y(:,i) = z;
    %             TT(i) = TT(i-1) + h;
    %             
    %             if ctr==0 && TT(i)>4.8
    %                 ctr = 1;
    %                 varStrict = z;
    %             end            
    %         else
    %             rejected = rejected + 1;
    %         end
    %         h = min(theta * (tol/err)^a * h,T-TT(i));
    %     end
    %     if mod(i,chunk)==0
    %         Y=[Y, Z];
    %         TT=[TT zeros(1,chunk)];
    %     end
    %  
    %     
    % end
    % Y=Y(:,1:i);
    % TT=TT(:,1:i);
    
    while TT(i)<T-5*eps
        err = tol + 1;
        while err>tol 
            [z,err] = variableRKMK45(vecField,action,Y(:,i),h);
            accepted = (err<tol);
            if accepted
                i = i + 1;
                Y(:,i) = z;
                TT(i) = TT(i-1) + h;

                if ctr==0
                    ctr = 1;
                    varStrict = z;
                end            
            else
                rejected = rejected + 1;
            end
            h = min(theta * (tol/err)^a * h,T-TT(i));
        end
        if mod(i,chunk)==0
            Y=[Y, Z];
            TT=[TT zeros(1,chunk)];
        end


    end
    Y=Y(:,1:i);
    TT=TT(:,1:i);


    Nsteps = i-1
    tt = linspace(0,T,Nsteps+1);
    dt = tt(2) - tt(1);
    zRK45 = zeros(length(z0),Nsteps);
    zRK45(:,1) = z0;
    count = 1;
    ts = 0;
    % for k=1:Nsteps
    %     zCF4 = FreeRK4SE3N(f,action,dt,zCF4);
    %     ts = ts + dt;
    %     if k*dt>4.8273 && ctr ==0
    %         ctr=1;
    %         zCF4Strict = zCF4;
    %     end
    % end

    for k=1:Nsteps
        zRK45(:,k+1) = RKMK5(vecField,action,zRK45(:,k),dt);
        ts = ts + dt;
    end


    % Solve with ODE45

    z0Ref = [q0;w0];
    odeFunc = @(t,z) [FuncQ(z);FuncW(z,L,m)];
    options = odeset('AbsTol',1e-12,'RelTol',1e-6);
    [timeSol,zC] = ode45(odeFunc, [0 T], z0Ref, options);
    zC = zC';
    zC = reorder(zC);
    zRef = zC(:,end);
    % 
    % hh = min(find((timeSol>4.8) == 1));
    % zCStrict = zC(:,hh);

    accVar(index) = norm(zRef-Y(:,end),2);
    accConst(index) = norm(zRef-zRK45(:,end),2);
    steps(index) = Nsteps;
    index = index + 1;
%     fprintf("Accuracy obtained for variable  stepsize: %10.2e\n",norm(zRef-Y(:,end),2))
%     fprintf("Accuracy obtained for constant stepsize: %10.2e\n",norm(zRef-zCF4,2));
%     fprintf("Number of rejected steps: %d2\n",rejected);
end
figure;
%subplot(1,2,1);
yyaxis left
semilogy(Prange,accConst,'r-o','linewidth',3,'Markersize',10);
hold on;
semilogy(Prange,accVar,'k-*','linewidth',3,'Markersize',10);
xlabel('Number of connected pendulums','FontSize',30);
ylabel('Accuracy at T=3','FontSize',30);
hold on;
ax = gca;
%grid on;
set(ax,'xcolor','k')
set(ax,'ycolor','k')
yyaxis right
ax = gca;
plot(Prange, steps,'b-s','linewidth',3,'Markersize',10);
set(ax,'ycolor','b')
ax.FontSize = 20;
%ax.YAxis.FontSize = 20;
ylabel('Number of time steps','FontSize',30);
%ax.YAxis.FontSize = 20;
%ax.XAxis.FontSize = 20;
%ax.YAxis.FontSize = 20;
% xlabel('Number of connected pendulums','FontSize',30);
% ylabel('Accuracy at T=3','FontSize',30)

legend('Accuracy RKMK5','Accuracy RKMK(5,4)','Number of time steps','FontSize',30);
%set(gca,'XTick',Prange,'XGrid','on')

% title("Accuracy against the number of pendulums",'FontSize',25);
%title("Accuracy against the number of pendulums",'FontSize',18);

%grid off


%% Comparison with step adaptation in ODE45
figure;
z0Ref = [q0;w0];
options = odeset('AbsTol',tol);
[timeSol,zC] = ode45(odeFunc, [0 T], z0Ref, options);
[err,times,sol] = variableStepComparison(vecField,action,z0,T,tol) ;
zC = zC';
zC=reorder(zC);
%subplot(1,2,2);
plot(TT(1:end-1),diff(TT),'k--','linewidth',2);
hold on;
plot(timeSol(1:end-1),diff(timeSol),'b:','linewidth',2);
hold on;
plot(timeSol(1:end-1),0*timeSol(1:end-1) + dt,'r-','linewidth',2);

legend('RKMK(5,4)','ODE45','RKMK5','FontSize',30,'Location','northwest');
xlim([0 T])
xlabel("Time",'fontsize',30);
ylabel("Stepsize",'fontsize',30);
ax = gca;
ax.XAxis.FontSize = 30;
ax.YAxis.FontSize = 30;
stringa = 'Comparison of stepsize variation with N='+string(P)+' pendulums';
%title(stringa,'FontSize',18);


% Prange = 20;
% for P = Prange
%     
%     L = rand(P,1)+0.5; %Set random lengths
%     m = rand(P,1)+0.5; %Set random masses
%     L = 0*L + Lref/P;%1; 
%     m = 0*m + 1;
% 
% %     [q0,w0,z0] = initializeSE3N(P);
%     [q0,w0,z0] = initializeSE3N_largeVariation(P);
% %     [q0,w0,z0] = initializeSE3N_smallVariation(P);
%     
%     Energy = @(q,w) 0.5*w'*assembleR(q,L,m)*w + potential(q,L,m);
% 
%     disp("Energy of this initial condition: "+num2str(Energy(q0,w0)));
% 
%     t0 = 0;
%     T = 3; 
%     N = 1000; 
%     time = linspace(t0,T,N); 
%     dt = time(2)-time(1);
% 
%     getq = @(v) extractq(v);
%     getw = @(v) extractw(v);
% 
%     f = @(v) fManiToAlgebra(getq(v),getw(v),L,m); 
%     action = @(B,input) actionSE3N(B,input); 
%     vecField = @(sigma,p) dexpinvSE3N(sigma,f(action(exponentialSE3N(sigma),p)));
% 
%     z = z0;
%     qC = zeros(3*P,N);
%     pC = qC;
% 
%     Len = zeros(3*P,1);
%     for i = 1:P
%         Len(3*i-2:3*i) = L(i)*ones(3,1);
%     end
%     Mat = diag(Len);
%     if P>1
%         for i = 3:3:3*(P-1)
%             Mat = Mat + diag(Len(1:3*P-i),-i);
%         end
%     end
%     qC(:,1) = q0;
%     pC(:,1) = Mat*q0;
% 
% 
%     %% COMPARISON BETWEEN CF34 AND CF4, i.e. VARIABLE STEPSIZE AGAINST CONSTANT ONE
% 
%     chunk=400;
%     Z=zeros(length(z0),chunk);
%     TT=zeros(1,chunk);
%     Y=Z;
%     Y(:,1)=z0;
% 
%     a = 1/4;
%     theta = 0.85;
%     i = 1;
%     h = T/100;
%     rejected = 0;
%     ctr = 0;
%     % while TT(i)<T-5*eps
%     %     err = tol + 1;
%     %     while err>tol    
%     %         [z,err] = VariableStepFreeRK4SE3N(f,action,h,Y(:,i));
%     %         accepted = (err<tol);
%     %         if accepted
%     %             i = i + 1;
%     %             Y(:,i) = z;
%     %             TT(i) = TT(i-1) + h;
%     %             
%     %             if ctr==0 && TT(i)>4.8
%     %                 ctr = 1;
%     %                 varStrict = z;
%     %             end            
%     %         else
%     %             rejected = rejected + 1;
%     %         end
%     %         h = min(theta * (tol/err)^a * h,T-TT(i));
%     %     end
%     %     if mod(i,chunk)==0
%     %         Y=[Y, Z];
%     %         TT=[TT zeros(1,chunk)];
%     %     end
%     %  
%     %     
%     % end
%     % Y=Y(:,1:i);
%     % TT=TT(:,1:i);
%     
%     while TT(i)<T-5*eps
%         err = tol + 1;
%         while err>tol 
%             [z,err] = variableRKMK45(vecField,action,Y(:,i),h);
%             accepted = (err<tol);
%             if accepted
%                 i = i + 1;
%                 Y(:,i) = z;
%                 TT(i) = TT(i-1) + h;
% 
%                 if ctr==0
%                     ctr = 1;
%                     varStrict = z;
%                 end            
%             else
%                 rejected = rejected + 1;
%             end
%             h = min(theta * (tol/err)^a * h,T-TT(i));
%         end
%         if mod(i,chunk)==0
%             Y=[Y, Z];
%             TT=[TT zeros(1,chunk)];
%         end
% 
% 
%     end
%     Y=Y(:,1:i);
%     TT=TT(:,1:i);
% 
% 
%     Nsteps = i-1
%     tt = linspace(0,T,Nsteps+1);
%     dt = tt(2) - tt(1);
%     zRK45 = zeros(length(z0),Nsteps);
%     zRK45(:,1) = z0;
%     count = 1;
%     ts = 0;
%     % for k=1:Nsteps
%     %     zCF4 = FreeRK4SE3N(f,action,dt,zCF4);
%     %     ts = ts + dt;
%     %     if k*dt>4.8273 && ctr ==0
%     %         ctr=1;
%     %         zCF4Strict = zCF4;
%     %     end
%     % end
% 
%     for k=1:Nsteps
%         zRK45(:,k+1) = RKMK5(vecField,action,zRK45(:,k),dt);
%         ts = ts + dt;
%     end
% 
% 
%     % Solve with ODE45
% 
%     z0Ref = [q0;w0];
%     odeFunc = @(t,z) [FuncQ(z);FuncW(z,L,m)];
%     options = odeset('AbsTol',1e-12,'RelTol',1e-6);
%     [timeSol,zC] = ode45(odeFunc, [0 T], z0Ref, options);
%     zC = zC';
%     zC = reorder(zC);
%     zRef = zC(:,end);
%     % 
%     % hh = min(find((timeSol>4.8) == 1));
%     % zCStrict = zC(:,hh);
% 
%     accVar(index) = norm(zRef-Y(:,end),2);
%     accConst(index) = norm(zRef-zRK45(:,end),2);
%     index = index + 1;
% 
% %     fprintf("Accuracy obtained for variable  stepsize: %10.2e\n",norm(zRef-Y(:,end),2))
% %     fprintf("Accuracy obtained for constant stepsize: %10.2e\n",norm(zRef-zCF4,2));
% %     fprintf("Number of rejected steps: %d2\n",rejected);
% end
% 
% figure;
% subplot(3,2,1);
% plot(TT,Y(1,:),'m-','linewidth',2);
% hold on;
% plot(TT,Y(2,:),'g-','linewidth',2);
% hold on;
% plot(TT,Y(3,:),'c-','linewidth',2);
% hold on;
% plot(TT,Y(4,:),'y-','linewidth',2);
% hold on;
% plot(TT,Y(5,:),'k-','linewidth',2);
% hold on;
% plot(TT,Y(6,:),'r-','linewidth',2);
% h = title("RKMK(5,4)");
% h.FontSize=30;
% ax = gca;
% ax.XAxis.FontSize = 40;
% ax.YAxis.FontSize = 40;
% 
% hSub = subplot(3,2,2); plot(1, nan,'m-', 1, nan,'g-', 1, nan,'c-', 1, nan,'y-',1, nan,'k-',1, nan,'r-', 'LineWidth', 3); set(hSub, 'Visible', 'off');
% 
% subplot(3,2,3);
% plot(timeSol,zC(1,:),'m-','linewidth',2);
% hold on;
% plot(timeSol,zC(2,:),'g-','linewidth',2);
% hold on;
% plot(timeSol,zC(3,:),'c-','linewidth',2);
% hold on;
% plot(timeSol,zC(4,:),'y-','linewidth',2);
% hold on;
% plot(timeSol,zC(5,:),'k-','linewidth',2);
% hold on;
% plot(timeSol,zC(6,:),'r-','linewidth',2);
% h = title("ODE45");
% h.FontSize=30;
% % legend('$q_1^{(1)}$','$q_1^{(2)}$','$q_1^{(3)}$','$\omega_1^{(1)}$','$\omega_1^{(2)}$',...
% %     '$\omega_1^{(3)}$','interpreter','latex','Location','northwest','FontSize',40)
% ax = gca;
% ax.XAxis.FontSize = 40;
% ax.YAxis.FontSize = 40;
% 
% subplot(3,2,5);
% plot(tt,zRK45(1,:),'m-','linewidth',2);
% hold on;
% plot(tt,zRK45(2,:),'g-','linewidth',2);
% hold on;
% plot(tt,zRK45(3,:),'c-','linewidth',2);
% hold on;
% plot(tt,zRK45(4,:),'y-','linewidth',2);
% hold on;
% plot(tt,zRK45(5,:),'k-','linewidth',2);
% hold on;
% plot(tt,zRK45(6,:),'r-','linewidth',2);
% xlabel("Time",'Fontsize',30);
% 
% h=title("RKMK5");
% h.FontSize=30;
% % legend('$q_1^{(1)}$','$q_1^{(2)}$','$q_1^{(3)}$','$\omega_1^{(1)}$','$\omega_1^{(2)}$',...
% %     '$\omega_1^{(3)}$','interpreter','latex','Location','northwest','FontSize',40)
% ax = gca;
% ax.XAxis.FontSize = 40;
% ax.YAxis.FontSize = 40;
% 
% legend(hSub,'$q_1^{(1)}$','$q_1^{(2)}$','$q_1^{(3)}$','$\omega_1^{(1)}$','$\omega_1^{(2)}$',...
%     '$\omega_1^{(3)}$','interpreter','latex','Location','northwest','FontSize',40)
% 
% figure;
% 
% subplot(3,2,1);
% plot(TT,Y(7,:),'m-','linewidth',2);
% hold on;
% plot(TT,Y(8,:),'g-','linewidth',2);
% hold on;
% plot(TT,Y(9,:),'c-','linewidth',2);
% hold on;
% plot(TT,Y(10,:),'y-','linewidth',2);
% hold on;
% plot(TT,Y(11,:),'k-','linewidth',2);
% hold on;
% plot(TT,Y(12,:),'r-','linewidth',2);
% h = title("RKMK(5,4)");
% h.FontSize=30;
% % legend('$q_2^{(1)}$','$q_2^{(2)}$','$q_2^{(3)}$','$\omega_2^{(1)}$','$\omega_2^{(2)}$',...
% %     '$\omega_2^{(3)}$','interpreter','latex','Location','northwest','FontSize',40)
% ax = gca;
% ax.XAxis.FontSize = 40;
% ax.YAxis.FontSize = 40;
% 
% hSub = subplot(3,2,2); plot(1, nan,'m-', 1, nan,'g-', 1, nan,'c-', 1, nan,'y-',1, nan,'k-',1, nan,'r-', 'LineWidth', 3); set(hSub, 'Visible', 'off');
% 
% subplot(3,2,3);
% plot(timeSol,zC(7,:),'m-','linewidth',2);
% hold on;
% plot(timeSol,zC(8,:),'g-','linewidth',2);
% hold on;
% plot(timeSol,zC(9,:),'c-','linewidth',2);
% hold on;
% plot(timeSol,zC(10,:),'y-','linewidth',2);
% hold on;
% plot(timeSol,zC(11,:),'k-','linewidth',2);
% hold on;
% plot(timeSol,zC(12,:),'r-','linewidth',2);
% h = title("ODE45");
% h.FontSize=30;
% % legend('$q_2^{(1)}$','$q_2^{(2)}$','$q_2^{(3)}$','$\omega_2^{(1)}$','$\omega_2^{(2)}$',...
% %     '$\omega_2^{(3)}$','interpreter','latex','Location','northwest','FontSize',40)
% ax = gca;
% ax.XAxis.FontSize = 40;
% ax.YAxis.FontSize = 40;
% 
% subplot(3,2,5);
% plot(tt,zRK45(7,:),'m-','linewidth',2);
% hold on;
% plot(tt,zRK45(8,:),'g-','linewidth',2);
% hold on;
% plot(tt,zRK45(9,:),'c-','linewidth',2);
% hold on;
% plot(tt,zRK45(10,:),'y-','linewidth',2);
% hold on;
% plot(tt,zRK45(11,:),'k-','linewidth',2);
% hold on;
% plot(tt,zRK45(12,:),'r-','linewidth',2);
% h = title("RKMK5");
% h.FontSize=30;
% xlabel("Time",'Fontsize',30);
% 
% ax = gca;
% ax.XAxis.FontSize = 40;
% ax.YAxis.FontSize = 40;
% 
% legend(hSub,'$q_2^{(1)}$','$q_2^{(2)}$','$q_2^{(3)}$','$\omega_2^{(1)}$','$\omega_2^{(2)}$',...
%     '$\omega_2^{(3)}$','interpreter','latex','Location','northwest','FontSize',40)
% 
%  
% 
% 
% %% zoomFinal
% figure;
% subplot(2,1,1);
% plot(tt,zRK45(11,:),'linewidth',2);
% hold on;
% plot(timeSol,zC(11,:),'linewidth',2);
% hold on;
% plot(TT,Y(11,:),'linewidth',2);
% hold on;
% xlim([2.1 T])
% xlabel("Time",'fontsize',40);
% % ss = ylabel("$\omega_1^{(2)}$",'interpreter','latex');
% ss = ylabel("\omega_1^{(2)}");
% ss.FontSize = 40;
% legend('RKMK5','ODE45','RKMK(5,4)','Location','northeast','FontSize',40);
% % h= title('Comparison of the component $\omega_2^{(2)}$','interpreter','latex','FontSize',40);
% h= title('Comparison of the component \omega_2^{(2)}');
% h.FontSize=40;
% ax = gca;
% ax.XAxis.FontSize = 40;
% ax.YAxis.FontSize = 40;
% subplot(2,1,2);
% plot(tt,zRK45(5,:),'linewidth',2);
% hold on;
% plot(timeSol,zC(5,:),'linewidth',2);
% hold on;
% plot(TT,Y(5,:),'linewidth',2);
% hold on;
% xlim([2.1 T])
% xlabel("Time",'fontsize',40);
% % ss = ylabel("$\omega_2^{(2)}$",'interpreter','latex');
% ss = ylabel("\omega_2^{(2)}");
% ss.FontSize = 40;
% legend('RKMK5','ODE45','RKMK(5,4)','Location','southeast','FontSize',40);
% % h = title('Comparison of the component $\omega_1^{(2)}$','interpreter','latex');
% h = title('Comparison of the component \omega_1^{(2)}');
% h.FontSize=40;
% ax = gca;
% ax.XAxis.FontSize = 40;
% ax.YAxis.FontSize = 40;
% 
% 
% 
% 
% %% stepAdaptP2 
% z0Ref = [q0;w0];
% options = odeset('AbsTol',tol);
% [timeSol,zC] = ode45(odeFunc, [0 T], z0Ref, options);
% [err,times,sol] = variableStepComparison(vecField,action,z0,T,tol) ;
% zC = zC';
% zC=reorder(zC);
% % subplot(1,2,2);
% plot(TT(1:end-1),diff(TT),'r-','linewidth',2);
% hold on;
% plot(timeSol(1:end-1),diff(timeSol),'k-','linewidth',2);
% hold on;
% plot(timeSol(1:end-1),0*timeSol(1:end-1) + dt,'m-','linewidth',2);
% 
% legend('RKMK(5,4)','ODE45','RKMK5','FontSize',40,'Location','northwest');
% xlim([0 T])
% xlabel("Time",'fontsize',40);
% ylabel("Stepsize",'fontsize',40);
% ax = gca;
% ax.XAxis.FontSize = 40;
% ax.YAxis.FontSize = 40;
% title('Comparison of stepsize variation with N=2 pendulums','FontSize',40);
% 
% 
% 
% 
% 
% % 
% % 
% % %% diffSol
% % z0Ref = [q0;w0];
% % options = odeset('AbsTol',tol);
% % [timeSol,zC] = ode45(odeFunc, [0 T], z0Ref, options);
% % [err,times,sol] = variableStepComparison(vecField,action,z0,T,tol) ;
% % zC = zC';
% % vecFieldw1 = zeros(1,length(timeSol));
% % vecFieldw2 = zeros(1,length(timeSol));
% % 
% % odeFunc = @(z) [FuncQ(z);FuncW(z,L,m)];
% % 
% % for i = 1:length(timeSol)
% %     vv = reorder(odeFunc(zC(:,i)));
% %     vecFieldw1(i) = vv(5);
% %     vecFieldw2(i) = vv(11);
% % end
% % 
% % figure;
% % plot(timeSol,vecFieldw1,'r-',timeSol,vecFieldw2,'k-','linewidth',2)
% % title("Two components of the vector field against time",'Fontsize',40)
% % xlabel("Time","Fontsize",40);
% % h = ylabel("$\dot{\omega}_i^{(2)}(t)$","interpreter","latex","Fontsize",40);
% % h.FontSize = 40;
% % 
% % l = legend('$\dot{\omega}_1^{(2)}(t)$','$\dot{\omega}_2^{(2)}(t)$','interpreter','latex',"Location","northwest");
% % l.FontSize = 60;
% % ylim([-110 110])
% % ax = gca;
% % ax.XAxis.FontSize = 40;
% % ax.YAxis.FontSize = 40;
