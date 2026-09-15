%% fits for the PK model of each anti-LSC drug

load('ProA_PK_data.mat');
load('ProA_PK_parameter_uncertainty.mat');
load('Digoxin_PK_data.mat');
load('Digoxin_PK_parameter_uncertainty.mat');
load('Oubain_PK_data.mat');
load('Ouabain_PK_parameter_uncertainty.mat');
load('Budesonide_PK_data.mat');
load('Budesonide_PK_parameter_uncertainty.mat');
load('Mometasone_PK_data.mat');
load('Mometasone_PK_parameter_uncertainty.mat')

ProA.units = (1/530.64968)*10^9; %mg/mL to nM
Dig.units = (1/780.938)*10^9;    %mg/mL to nM
Oua.units = (1/584.6525)*10^9;   %mg/mL to nM
Bud.units = (1/430.534)*10^9;    %mg/mL to nM
Mom.units = (1/521.4)*10^9;      %mg/mL to nM

%set up data
ProA.PK_dataTime = ProA_PK_data(:,2); %in min
ProA.PK_dataConc = ProA_PK_data(:,3); %in mg/mL
ProA.V_popSA = ProA_PK_parameter_uncertainty(:,1);
ProA.k_popSA = ProA_PK_parameter_uncertainty(:,2);
ProA.k12_popSA = ProA_PK_parameter_uncertainty(:,3);
ProA.k21_popSA = ProA_PK_parameter_uncertainty(:,4);

Dig.PK_dataTime = Digoxin_PK_data(:,2); %in min
Dig.PK_dataConc = Digoxin_PK_data(:,3); %in mg/mL
Dig.V_popSA = Digoxin_PK_parameter_uncertainty(:,1);
Dig.k_popSA = Digoxin_PK_parameter_uncertainty(:,2);
Dig.k12_popSA = Digoxin_PK_parameter_uncertainty(:,3);
Dig.k21_popSA = Digoxin_PK_parameter_uncertainty(:,4);

Oua.PK_dataTime = Oubain_PK_data(:,2); %in min
Oua.PK_dataConc = Oubain_PK_data(:,3); %in mg/mL
Oua.V_popSA = Ouabain_PK_parameter_uncertainty(:,1);
Oua.k_popSA = Ouabain_PK_parameter_uncertainty(:,2);
Oua.k12_popSA = Ouabain_PK_parameter_uncertainty(:,3);
Oua.k21_popSA = Ouabain_PK_parameter_uncertainty(:,4);
Oua.k13_popSA = Ouabain_PK_parameter_uncertainty(:,5);
Oua.k31_popSA = Ouabain_PK_parameter_uncertainty(:,6);

Bud.PK_dataTime = Budesonide_PK_data(:,2); %in min
Bud.PK_dataConc = Budesonide_PK_data(:,3); %in mg/mL
Bud.V_popSA = Budesonide_PK_parameter_uncertainty(:,1);
Bud.k_popSA = Budesonide_PK_parameter_uncertainty(:,2);
Bud.k12_popSA = Budesonide_PK_parameter_uncertainty(:,3);
Bud.k21_popSA = Budesonide_PK_parameter_uncertainty(:,4);

Mom.PK_dataTime = Mometasone_PK_data(:,2); %in min
Mom.PK_dataConc = Mometasone_PK_data(:,3); %in mg/mL
Mom.V_popSA = Mometasone_PK_parameter_uncertainty(:,1);
Mom.k_popSA = Mometasone_PK_parameter_uncertainty(:,2);
Mom.k12_popSA = Mometasone_PK_parameter_uncertainty(:,3);
Mom.k21_popSA = Mometasone_PK_parameter_uncertainty(:,4);
Mom.k13_popSA = Mometasone_PK_parameter_uncertainty(:,5);
Mom.k31_popSA = Mometasone_PK_parameter_uncertainty(:,6);

%parameters from PK fit
ProA.D = 1.2;             %in mg
ProA.t_inf = 70;          %in min
ProA.compartmentNb = 2;
ProA.V_pop = 24.42*1000;  %in mL
ProA.k_pop = 0.78/60;     %in 1/min
ProA.k12_pop = 1.02/60;   %in 1/min
ProA.k21_pop = 0.0402/60; %in 1/min

Dig.D = 1;                %in mg
Dig.t_inf = 10;           %in min
Dig.compartmentNb = 2;
Dig.V_pop = 38.41*1000;   %in mL
Dig.k_pop = 0.156/60;     %in 1/min
Dig.k12_pop = 0.66/60;    %in 1/min
Dig.k21_pop = 0.096/60;   %in 1/min
Dig.omega_k12 = 0.38;
Dig.omega_k21 = 0.59;

Oua.D = 0.5;              %in mg
Oua.compartmentNb = 3;
Oua.V_pop = 41.93*1000;   %in mL
Oua.k_pop = 0.96/60;      %in 1/min
Oua.k12_pop = 5.10/60;    %in 1/min
Oua.k21_pop = 1.50/60;    %in 1/min
Oua.k13_pop = 2.16/60;    %in 1/min
Oua.k31_pop = 0.108/60;   %in 1/min

Bud.D = 0.5;              %in mg
Bud.t_inf = 10;           %in min
Bud.compartmentNb = 2;
Bud.V_pop = 33.18*1000;   %in mL
Bud.k_pop = 2.34/60;      %in 1/min
Bud.k12_pop = 4.68/60;    %in 1/min
Bud.k21_pop = 1.14/60;    %in 1/min

Mom.D = 0.4;              %in mg
Mom.t_inf = 1;            %in min
Mom.compartmentNb = 3;
Mom.V_pop = 45.58*1000;   %in mL
Mom.k_pop = 0.84/60;      %in 1/min
Mom.k12_pop = 0.09/60;    %in 1/min
Mom.k21_pop = 0.0444/60;  %in 1/min
Mom.k13_pop = 0.432/60;   %in 1/min
Mom.k31_pop = 0.444/60;   %in 1/min

%% Simulating PK model

%Proscillaridin A
ProA.PK_simTime = 0:1:4320; %in min
ProA.In = zeros(1,size(ProA.PK_simTime,2));
for i = 1:size(ProA.In,2)
    if ProA.PK_simTime(1,i) <= ProA.t_inf 
        ProA.In(1,i) = ProA.D/ProA.t_inf; 
    else
        ProA.In(1,i) = 0;
    end
end
ProA.PK_initialConds = [0;0];

sol = ode45(@(t,A) fitting_ProA_PK(t,A,ProA), ProA.PK_simTime, ProA.PK_initialConds);
ProA.Ac = deval(sol,ProA.PK_simTime,1);
ProA.Cc = ProA.Ac/ProA.V_pop;

%Digoxin
Dig.PK_simTime = 0:1:2880; %in min
Dig.In = zeros(1,size(Dig.PK_simTime,2));
for i = 1:size(Dig.In,2)
    if Dig.PK_simTime(1,i) <= Dig.t_inf
        Dig.In(1,i) = Dig.D/Dig.t_inf;
    else
        Dig.In(1,i) = 0;
    end
end
Dig.PK_initialConds = [0 0];

sol = ode45(@(t,A) fitting_Digoxin_PK(t,A,Dig), Dig.PK_simTime, Dig.PK_initialConds);
Dig.Ac = deval(sol,Dig.PK_simTime,1);
Dig.Cc = Dig.Ac/Dig.V_pop;

%Ouabain
Oua.PK_simTime = 0:1:2880; %in min
Oua.PK_initialConds = [Oua.D;0;0];

sol = ode45(@(t,A) fitting_Ouabain_PK(t,A,Oua), Oua.PK_simTime, Oua.PK_initialConds);
Oua.Ac = deval(sol,Oua.PK_simTime,1);
Oua.Cc = Oua.Ac/Oua.V_pop;

%Budesonide
Bud.PK_simTime = 0:1:2880; %in min
Bud.In = zeros(1,size(Bud.PK_simTime,2));
for i = 1:size(Bud.In,2)
    if Bud.PK_simTime(1,i) <= Bud.t_inf
        Bud.In(1,i) = Bud.D/Bud.t_inf;
    else
        Bud.In(1,i) = 0;
    end
end
Bud.PK_initialConds = [0;0];

sol = ode45(@(t,A) fitting_Budesonide_PK(t,A,Bud), Bud.PK_simTime, Bud.PK_initialConds);
Bud.Ac = deval(sol,Bud.PK_simTime,1);
Bud.Cc = Bud.Ac/Bud.V_pop; %Plasma concentration

%Mometasone
Mom.PK_simTime = 0:1:2880; %in min
Mom.PK_initialConds = [Mom.D;0;0];

sol = ode45(@(t,A) fitting_Mometasone_PK(t,A,Mom), Mom.PK_simTime, Mom.PK_initialConds);
Mom.Ac = deval(sol,Mom.PK_simTime,1);
Mom.Cc = Mom.Ac/Mom.V_pop;

%% 95% confidence bands of PK fit

%Proscillaridin A
ProA.Ac_SA = zeros(size(ProA_PK_parameter_uncertainty,1),size(ProA.PK_simTime,2));

for i = 1: size(ProA_PK_parameter_uncertainty,1)
    ProA.V_pop = ProA.V_popSA(i);     %in mL
    ProA.k_pop = ProA.k_popSA(i);     %in 1/min
    ProA.k12_pop = ProA.k12_popSA(i); %in 1/min
    ProA.k21_pop = ProA.k21_popSA(i); %in 1/min    
    
    sol = ode45(@(t,A) fitting_ProA_PK(t,A,ProA), ProA.PK_simTime, ProA.PK_initialConds);
    ProA.Ac_SA(i,:) = deval(sol,ProA.PK_simTime,1);
    ProA.Cc_SA(i,:) = ProA.Ac_SA(i,:)/ProA.V_pop;
end

ProA.PK_ci = zeros(2,size(ProA.Cc_SA,2));
for i = 1:size(ProA.Cc_SA,2)
   ProA.PK_ci(:,i) = prctile(ProA.Cc_SA(:,i),[2.5,97.5]);
end

%Digoxin
Dig.Ac_SA = zeros(size(Digoxin_PK_parameter_uncertainty,1),size(Dig.PK_simTime,2));

for i = 1: size(Digoxin_PK_parameter_uncertainty,1)
    Dig.V_pop = Dig.V_popSA(i);     %in mL
    Dig.k_pop = Dig.k_popSA(i);     %in 1/min
    Dig.k12_pop = Dig.k12_popSA(i); %in 1/min
    Dig.k21_pop = Dig.k21_popSA(i); %in 1/min    
    
    sol = ode45(@(t,A) fitting_Digoxin_PK(t,A,Dig), Dig.PK_simTime, Dig.PK_initialConds);
    Dig.Ac_SA(i,:) = deval(sol,Dig.PK_simTime,1);
    Dig.Cc_SA(i,:) = Dig.Ac_SA(i,:)/Dig.V_pop;
end

Dig.PK_ci = zeros(2,size(Dig.Cc_SA,2));
for i = 1:size(Dig.Cc_SA,2)
   Dig.PK_ci(:,i) = prctile(Dig.Cc_SA(:,i),[2.5,97.5]);
end

%Ouabain
Oua.Ac_SA = zeros(size(Ouabain_PK_parameter_uncertainty,1),size(Oua.PK_simTime,2));

for i = 1: size(Ouabain_PK_parameter_uncertainty,1)
    Oua.V_pop = Oua.V_popSA(i);     %in mL
    Oua.k_pop = Oua.k_popSA(i);     %in 1/min
    Oua.k12_pop = Oua.k12_popSA(i); %in 1/min
    Oua.k21_pop = Oua.k21_popSA(i); %in 1/min  
    Oua.k13_pop = Oua.k13_popSA(i); %in 1/min
    Oua.k31_pop = Oua.k31_popSA(i); %in 1/min 
    
    sol = ode45(@(t,A) fitting_Ouabain_PK(t,A,Oua), Oua.PK_simTime, Oua.PK_initialConds);
    Oua.Ac_SA(i,:) = deval(sol,Oua.PK_simTime,1);
    Oua.Cc_SA(i,:) = Oua.Ac_SA(i,:)/Oua.V_pop; 
end

Oua.PK_ci = zeros(2,size(Oua.Cc_SA,2));
for i = 1:size(Oua.Cc_SA,2)
   Oua.PK_ci(:,i) = prctile(Oua.Cc_SA(:,i),[2.5,97.5]);
end

%Budesonide
Bud.Ac_SA = zeros(size(Budesonide_PK_parameter_uncertainty,1),size(Bud.PK_simTime,2));

for i = 1: size(Budesonide_PK_parameter_uncertainty,1)
    Bud.V_pop = Bud.V_popSA(i);     %in mL
    Bud.k_pop = Bud.k_popSA(i);     %in 1/min
    Bud.k12_pop = Bud.k12_popSA(i); %in 1/min
    Bud.k21_pop = Bud.k21_popSA(i); %in 1/min  
    
    sol = ode45(@(t,A) fitting_Budesonide_PK(t,A,Bud), Bud.PK_simTime, Bud.PK_initialConds);
    Bud.Ac_SA(i,:) = deval(sol,Bud.PK_simTime,1);
    Bud.Cc_SA(i,:) = Bud.Ac_SA(i,:)/Bud.V_pop;
end

Bud.PK_ci = zeros(2,size(Bud.Cc_SA,2));
for i = 1:size(Bud.Cc_SA,2)
   Bud.PK_ci(:,i) = prctile(Bud.Cc_SA(:,i),[2.5,97.5]);
end

%Mometasone
Mom.Ac_SA = zeros(size(Mometasone_PK_parameter_uncertainty,1),size(Mom.PK_simTime,2));

for i = 1: size(Mometasone_PK_parameter_uncertainty,1)
    Mom.V_pop = Mom.V_popSA(i);     %in mL
    Mom.k_pop = Mom.k_popSA(i);     %in 1/min
    Mom.k12_pop = Mom.k12_popSA(i); %in 1/min
    Mom.k21_pop = Mom.k21_popSA(i); %in 1/min  
    Mom.k13_pop = Mom.k13_popSA(i); %in 1/min
    Mom.k31_pop = Mom.k31_popSA(i); %in 1/min 
    
    sol = ode45(@(t,A) fitting_Mometasone_PK(t,A,Mom), Mom.PK_simTime, Mom.PK_initialConds);
    Mom.Ac_SA(i,:) = deval(sol,Mom.PK_simTime,1);
    Mom.Cc_SA(i,:) = Mom.Ac_SA(i,:)/Mom.V_pop;
end

Mom.PK_ci = zeros(2,size(Mom.Cc_SA,2));
for i = 1:size(Mom.Cc_SA,2)
   Mom.PK_ci(:,i) = prctile(Mom.Cc_SA(:,i),[2.5,97.5]);
end

%% Save work

% save('ProA_PK.mat', '-struct', 'ProA');
% save('Digoxin_PK.mat', '-struct', 'Dig');
% save('Ouabain_PK.mat', '-struct', 'Oua');
% save('Budesonide_PK.mat', '-struct', 'Bud');
% save('Mometasone_PK.mat', '-struct', 'Mom');

%% Figure S1: Population pharmacokinetics for each candidate drug

figure
tiledlayout(2,3,'TileSpacing','compact');

%Proscillaridin A
nexttile
hold on 
scatter(ProA.PK_dataTime/60, ProA.PK_dataConc.*ProA.units,'MarkerEdgeColor','#001253','MarkerFaceColor','#001253');
plot(ProA.PK_simTime/60,ProA.Cc.*ProA.units,'Color','#001253','LineWidth',1.5); % plotting curve in nM per hour
patch([ProA.PK_simTime/60, fliplr(ProA.PK_simTime/60)],[ProA.PK_ci(1,:)*ProA.units, fliplr(ProA.PK_ci(2,:)*ProA.units)],1,'facecolor','#001253', 'edgecolor', 'none', 'facealpha', 0.1); %CI
hold off
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
set(gca,'xtick',0:6:24);
xlabel('Time (hours)')
ylabel('Concentration (nM)')
title('Proscilaridin A','FontSize',20)
xlim([0,24])

%Digoxin
nexttile
hold on 
scatter(Dig.PK_dataTime/60, Dig.PK_dataConc.*Dig.units,'MarkerEdgeColor','#001253','MarkerFaceColor','#001253');
plot(Dig.PK_simTime/60,Dig.Cc.*Dig.units,'Color','#001253','LineWidth',1.5); % plotting curve in nM per hour
patch([Dig.PK_simTime/60, fliplr(Dig.PK_simTime/60)],[Dig.PK_ci(1,:)*Dig.units, fliplr(Dig.PK_ci(2,:)*Dig.units)],1,'facecolor','#001253', 'edgecolor', 'none', 'facealpha', 0.1); %CI
hold off
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
set(gca,'xtick',0:6:24);
xlabel('Time (hours)')
ylabel('Concentration (nM)')
title('Digoxin','FontSize',20)
xlim([0,24])

%Ouabain
nexttile
hold on 
scatter(Oua.PK_dataTime/60, Oua.PK_dataConc.*Oua.units,'MarkerEdgeColor','#001253','MarkerFaceColor','#001253');
plot(Oua.PK_simTime/60,Oua.Cc.*Oua.units,'Color','#001253','LineWidth',1.5); % plotting curve
patch([Oua.PK_simTime/60, fliplr(Oua.PK_simTime/60)],[Oua.PK_ci(1,:)*Oua.units, fliplr(Oua.PK_ci(2,:)*Oua.units)],1,'facecolor','#001253', 'edgecolor', 'none', 'facealpha', 0.1); %CI
hold off
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
set(gca,'xtick',0:6:24,'yscale','log');
xlabel('Time (hours)')
ylabel('Concentration (nM)')
title('Ouabain','FontSize',20)
xlim([0,24])

%Budesonide
nexttile
hold on 
scatter(Bud.PK_dataTime/60, Bud.PK_dataConc.*Bud.units,'MarkerEdgeColor','#001253','MarkerFaceColor','#001253');
plot(Bud.PK_simTime/60,Bud.Cc.*Bud.units,'Color','#001253','LineWidth',1.5); % plotting curve
patch([Bud.PK_simTime/60, fliplr(Bud.PK_simTime/60)],[Bud.PK_ci(1,:)*Bud.units, fliplr(Bud.PK_ci(2,:)*Bud.units)],1,'facecolor','#001253', 'edgecolor', 'none', 'facealpha', 0.1); %CI
hold off
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
set(gca,'xtick',0:6:24);
xlabel('Time (hours)')
ylabel('Concentration (nM)')
title('Budesonide','FontSize',20)
xlim([0,24])

%Mometasone
nexttile
hold on 
scatter(Mom.PK_dataTime/60, Mom.PK_dataConc.*Mom.units,'MarkerEdgeColor','#001253','MarkerFaceColor','#001253');
plot(Mom.PK_simTime/60,Mom.Cc.*Mom.units,'Color','#001253','LineWidth',1.5); % plotting curve
patch([Mom.PK_simTime/60, fliplr(Mom.PK_simTime/60)],[Mom.PK_ci(1,:)*Mom.units, fliplr(Mom.PK_ci(2,:)*Mom.units)],1,'facecolor','#001253', 'edgecolor', 'none', 'facealpha', 0.1); %CI
hold off
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
set(gca,'xtick',0:6:24);
xlabel('Time (hours)')
ylabel('Concentration (nM)')
title('Mometasone','FontSize',20)
xlim([0,24])

