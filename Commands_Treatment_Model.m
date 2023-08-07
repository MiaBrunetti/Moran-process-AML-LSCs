%% Commands for treatment model

%load work
ProA = load('ProA_PKPD.mat');
Dig = load('Digoxin_PKPD.mat');
Oua = load('Ouabain_PKPD.mat');
Bud = load('Budesonide_PKPD.mat');
Mom = load('Mometasone_PKPD.mat');

modelfn_Tox = @(beta,C) beta(1)*(C.^beta(3))./(beta(2).^beta(3) + C.^beta(3));
modelfn_Viab = @(beta,C) beta(2) - (beta(2)*C.^beta(4))./(beta(3).^beta(4) + C.^beta(4)) + beta(1);

%parameters of Moran process
sim_time = 365*24*60;                         %in min
Nd = 50000;                                   %total number of stem cells
r = 28*24*60;                                 %division rate of one stem cell in minutes: 1 HCS divides every r minutes (28 days)
div = round((1/(r/Nd))*100)/100;              %number of cell divisions in a minute
maxNoEvents = sim_time*div;                   %number of events/divisions in the Moran Process for the simulated time
timeofDiv = linspace(0,sim_time,maxNoEvents); %time at which each division occurs for a full year

%setting doses 
CarGly_doses_nM = [50 30 20 10]; %in nM
Glu_doses_nM = [25 10 1.5 0.25]; %in nM

%% Daily IV bolus treatment simulation

%Proscillaridin A
ProA.doses_nM = CarGly_doses_nM; % in nM
ProA.doses = (ProA.doses_nM/ProA.units)*ProA.V_pop; % in mg

ProA.NumberAdmins = 365;            
ProA.TimeFirstAdmin = 0;            %time of first administration
ProA.DosingInterval = 24*60;        %dosing interval in min
ProA.TotalTime =...
    [ProA.TimeFirstAdmin ProA.NumberAdmins*ProA.DosingInterval]; %in min
ProA.TreatmentTime =...
    ProA.TimeFirstAdmin:1:ProA.NumberAdmins*ProA.DosingInterval; %in min
ProA.AdministrationTimes =...
    linspace(ProA.TimeFirstAdmin,ProA.TimeFirstAdmin+(ProA.DosingInterval*(ProA.NumberAdmins-1)),ProA.NumberAdmins);

ProA.TreatmentAc = zeros(size(ProA.doses,2),size(timeofDiv,2));
ProA.TreatmentCc = zeros(size(ProA.TreatmentAc));
ProA.TreatmentToxicity = zeros(size(ProA.TreatmentAc));
ProA.Treatment_vHSC = zeros(size(ProA.TreatmentAc));
ProA.Treatment_vLSC = zeros(size(ProA.TreatmentAc));

for i = 1:size(ProA.doses,2)
    ProA.Dose = ProA.doses(1,i); %in mg
    ProA.initial_conds = [ProA.Dose 0];
    sol = simulation_PKPD_model(ProA); %solve ODE system to get treatment PK
    ProA.TreatmentAc(i,:) = deval(sol,timeofDiv,1);
    for j = 2:ProA.NumberAdmins
        idx = find(timeofDiv == ProA.AdministrationTimes(j));
        if ~isempty(idx)
            ProA.TreatmentAc(i,idx) = ProA.TreatmentAc(i,idx) + ProA.Dose/2;
        end
    end
    ProA.TreatmentCc(i,:) = (ProA.TreatmentAc(i,:)/ProA.V_pop)*ProA.units; %mg to nM
    ProA.TreatmentToxicity(i,:) = modelfn_Tox(ProA.Tox_Para_fit,ProA.TreatmentCc(i,:));
    ProA.Treatment_vHSC(i,:) = modelfn_Viab(ProA.vHSC_Para_fit,ProA.TreatmentCc(i,:)); 
    ProA.Treatment_vLSC(i,:) = modelfn_Viab(ProA.vLSC_Para_fit,ProA.TreatmentCc(i,:)); 
    
    disp(['ProA Loop #' num2str(i) ''])
end

%Digoxin
Dig.doses_nM = CarGly_doses_nM; %in nM
Dig.doses = (Dig.doses_nM/Dig.units)*Dig.V_pop; %in mg

Dig.NumberAdmins = 365;
Dig.TimeFirstAdmin = 0;            %time of first administration
Dig.DosingInterval = 24*60;        %dosing interval in min
Dig.TotalTime =...
    [Dig.TimeFirstAdmin Dig.NumberAdmins*Dig.DosingInterval]; %in min
Dig.TreatmentTime =...
    Dig.TimeFirstAdmin:1:Dig.NumberAdmins*Dig.DosingInterval; %in min
Dig.AdministrationTimes =...
    linspace(Dig.TimeFirstAdmin,Dig.TimeFirstAdmin+(Dig.DosingInterval*(Dig.NumberAdmins-1)),Dig.NumberAdmins);

Dig.TreatmentAc = zeros(size(Dig.doses,2),size(timeofDiv,2));
Dig.TreatmentCc = zeros(size(Dig.TreatmentAc));
Dig.TreatmentToxicity = zeros(size(Dig.TreatmentAc));
Dig.Treatment_vHSC = zeros(size(Dig.TreatmentAc));
Dig.Treatment_vLSC = zeros(size(Dig.TreatmentAc));

for i = 1:size(Dig.doses,2)
    Dig.Dose = Dig.doses(1,i); %in mg
    Dig.initial_conds = [Dig.Dose 0];
    
    sol = simulation_PKPD_model(Dig); %solve ODE system for treatment PK
    Dig.TreatmentAc(i,:) = deval(sol,timeofDiv,1);
    for j = 2:Dig.NumberAdmins
        idx = find(timeofDiv == Dig.AdministrationTimes(j));
        if ~isempty(idx)
            Dig.TreatmentAc(i,idx) = Dig.TreatmentAc(i,idx) + Dig.Dose/2;
        end
    end
    Dig.TreatmentCc(i,:) = (Dig.TreatmentAc(i,:)/Dig.V_pop)*Dig.units; %mg to nM
    Dig.TreatmentToxicity(i,:) = modelfn_Tox(Dig.Tox_Para_fit,Dig.TreatmentCc(i,:));
    Dig.Treatment_vHSC(i,:) = modelfn_Viab(Dig.vHSC_Para_fit,Dig.TreatmentCc(i,:)); 
    Dig.Treatment_vLSC(i,:) = modelfn_Viab(Dig.vLSC_Para_fit,Dig.TreatmentCc(i,:));
    
    disp(['Digoxin Loop #' num2str(i) ''])
end

%Ouabain
Oua.doses_nM = CarGly_doses_nM; %in nM
Oua.doses = (Oua.doses_nM/Oua.units)*Oua.V_pop; %in mg

Oua.NumberAdmins = 365;
Oua.TimeFirstAdmin = 0;            %time of first administration
Oua.DosingInterval = 24*60;        %dosing interval in min
Oua.TotalTime =...
    [Oua.TimeFirstAdmin Oua.NumberAdmins*Oua.DosingInterval]; %in min
Oua.TreatmentTime =...
    Oua.TimeFirstAdmin:1:Oua.NumberAdmins*Oua.DosingInterval; %in min
Oua.AdministrationTimes =...
    linspace(Oua.TimeFirstAdmin,Oua.TimeFirstAdmin+(Oua.DosingInterval*(Oua.NumberAdmins-1)),Oua.NumberAdmins);

Oua.TreatmentAc = zeros(size(Oua.doses,2),size(timeofDiv,2));
Oua.TreatmentCc = zeros(size(Oua.TreatmentAc));
Oua.TreatmentToxicity = zeros(size(Oua.TreatmentAc));
Oua.Treatment_vHSC = zeros(size(Oua.TreatmentAc));
Oua.Treatment_vLSC = zeros(size(Oua.TreatmentAc));

for i = 1:size(Oua.doses,2)   
    Oua.Dose = Oua.doses(1,i); %in mg
    Oua.initial_conds = [Oua.Dose 0 0];
    
    sol = simulation_PKPD_model(Oua); %solve ODE system to get treatment PK
    Oua.TreatmentAc(i,:) = deval(sol,timeofDiv,1);
    for j = 2:Oua.NumberAdmins
        [val,idx] = min(abs(timeofDiv - Oua.AdministrationTimes(j)));
        Oua.TreatmentAc(i,idx) = Oua.TreatmentAc(i,idx-1) + Oua.Dose;
    end
    
    Oua.TreatmentCc(i,:) = (Oua.TreatmentAc(i,:)/Oua.V_pop)*Oua.units; %mg to nM 
    Oua.TreatmentToxicity(i,:) = modelfn_Tox(Oua.Tox_Para_fit,Oua.TreatmentCc(i,:));
    Oua.Treatment_vHSC(i,:) = modelfn_Viab(Oua.vHSC_Para_fit,Oua.TreatmentCc(i,:));  
    Oua.Treatment_vLSC(i,:) = modelfn_Viab(Oua.vLSC_Para_fit,Oua.TreatmentCc(i,:));  
    
    disp(['Ouabain Loop #' num2str(i) ''])
end

%Budesonide
Bud.doses_nM = Glu_doses_nM; %in nM
Bud.doses = (Bud.doses_nM/Bud.units)*Bud.V_pop; %in mg

Bud.NumberAdmins = 365;
Bud.TimeFirstAdmin = 0;            %time of first administration
Bud.DosingInterval = 24*60;        %dosing interval in min
Bud.TotalTime =...
    [Bud.TimeFirstAdmin Bud.NumberAdmins*Bud.DosingInterval]; %in min
Bud.TreatmentTime =...
    Bud.TimeFirstAdmin:1:Bud.NumberAdmins*Bud.DosingInterval; %in min
Bud.AdministrationTimes =...
    linspace(Bud.TimeFirstAdmin,Bud.TimeFirstAdmin+(Bud.DosingInterval*(Bud.NumberAdmins-1)),Bud.NumberAdmins);

Bud.TreatmentAc = zeros(size(Bud.doses,2),size(timeofDiv,2));
Bud.TreatmentCc = zeros(size(Bud.TreatmentAc));
Bud.Treatment_vHSC = zeros(size(Bud.TreatmentAc));
Bud.Treatment_vLSC = zeros(size(Bud.TreatmentAc));

for i = 1:size(Bud.doses,2)   
    Bud.Dose = Bud.doses(1,i); %in mg
    Bud.initial_conds = [Bud.Dose 0];
    
    sol = simulation_PKPD_model(Bud); %solve ODE system to get treatment PK
    Bud.TreatmentAc(i,:) = deval(sol,timeofDiv,1);
    for j = 2:Bud.NumberAdmins
        idx = find(timeofDiv == Bud.AdministrationTimes(j));
        if ~isempty(idx)
            Bud.TreatmentAc(i,idx) = Bud.TreatmentAc(i,idx) + Bud.Dose/2;
        end
    end
    Bud.TreatmentCc(i,:) = (Bud.TreatmentAc(i,:)/Bud.V_pop)*Bud.units; % in mg to nM
    Bud.Treatment_vHSC(i,:) = modelfn_Viab(Bud.vHSC_Para_fit,Bud.TreatmentCc(i,:)); 
    Bud.Treatment_vLSC(i,:) = modelfn_Viab(Bud.vLSC_Para_fit,Bud.TreatmentCc(i,:)); 
    
    disp(['Budesonide Loop #' num2str(i) ''])
end

%Mometasone
Mom.doses_nM = Glu_doses_nM; %in nM
Mom.doses = (Mom.doses_nM/Mom.units)*Mom.V_pop; %in mg

Mom.NumberAdmins = 365;
Mom.TimeFirstAdmin = 0;            %time of first administration
Mom.DosingInterval = 24*60;        %in min
Mom.TotalTime =...
    [Mom.TimeFirstAdmin Mom.NumberAdmins*Mom.DosingInterval]; %in min
Mom.TreatmentTime =...
    Mom.TimeFirstAdmin:1:Mom.NumberAdmins*Mom.DosingInterval; %in min
Mom.AdministrationTimes =...
    linspace(Mom.TimeFirstAdmin,Mom.TimeFirstAdmin+(Mom.DosingInterval*(Mom.NumberAdmins-1)),Mom.NumberAdmins);

Mom.TreatmentAc = zeros(size(Mom.doses,2),size(timeofDiv,2));
Mom.TreatmentCc = zeros(size(Mom.TreatmentAc));
Mom.Treatment_vHSC = zeros(size(Mom.TreatmentAc));
Mom.Treatment_vLSC = zeros(size(Mom.TreatmentAc));

for i = 1:size(Mom.doses,2)
    Mom.Dose = Mom.doses(1,i); %in mg
    Mom.initial_conds = [Mom.Dose 0 0];
    
    sol = simulation_PKPD_model(Mom); %solve ODE system to get treatment PK 
    Mom.TreatmentAc(i,:) = deval(sol,timeofDiv,1);
    for j = 2:Mom.NumberAdmins
        idx = find(timeofDiv == Mom.AdministrationTimes(j));
        if ~isempty(idx)
            Mom.TreatmentAc(i,idx) = Mom.TreatmentAc(i,idx) + Mom.Dose/2;
        end
    end
    Mom.TreatmentCc(i,:) = (Mom.TreatmentAc(i,:)/Mom.V_pop)*Mom.units; 
    Mom.Treatment_vHSC(i,:) =  modelfn_Viab(Mom.vHSC_Para_fit,Mom.TreatmentCc(i,:)); 
    Mom.Treatment_vLSC(i,:) = modelfn_Viab(Mom.vLSC_Para_fit,Mom.TreatmentCc(i,:)); 
    
    disp(['Mometasone Loop #' num2str(i) ''])
end

%%  Finding steady state toxicity

d1 = 33927; %day 19 (20th adminitration)
d2 = 20*24*60*div; %day 20

%Proscillaridin A
ProA.CSS = zeros(length(CarGly_doses_nM),1);
ProA.ToxSS = zeros(length(CarGly_doses_nM),1);

for i = 1:length(CarGly_doses_nM)
    ProA.CSS(i,1) = trapz(ProA.TreatmentCc(i,d1:d2))/(60*24);
    ProA.ToxSS(i,1) = modelfn_Tox(ProA.Tox_Para_fit,ProA.CSS(i,1));
end

%Digoxin
Dig.CSS = zeros(length(CarGly_doses_nM),1);
Dig.ToxSS = zeros(length(CarGly_doses_nM),1);

for i = 1:length(CarGly_doses_nM)
    Dig.CSS(i,1) = trapz(Dig.TreatmentCc(i,d1:d2))/(60*24);
    Dig.ToxSS(i,1) = modelfn_Tox(Dig.Tox_Para_fit,Dig.CSS(i,1));
end

%Ouabain
Oua.CSS = zeros(length(CarGly_doses_nM),1);
Oua.ToxSS = zeros(length(CarGly_doses_nM),1);

for i = 1:length(CarGly_doses_nM)
    Oua.CSS(i,1) = trapz(Oua.TreatmentCc(i,d1:d2))/(60*24);
    Oua.ToxSS(i,1) = modelfn_Tox(Oua.Tox_Para_fit,Oua.CSS(i,1));
end

%% Save work

% save('ProA_TreatmentModel.mat', '-struct', 'ProA');
% save('Digoxin_TreatmentModel.mat', '-struct', 'Dig');
% save('Ouabain_TreatmentModel.mat', '-struct', 'Oua');
% save('Budesonide_TreatmentModel.mat', '-struct', 'Bud');
% save('Mometasone_TreatmentModel.mat', '-struct', 'Mom');

%% Figure 4A: Predicted PKPD responses of candidate cardiac glycosides

figure
tiledlayout(2,3,'TileSpacing','compact');
c = char('#702A8C','#BF2669','#FF7326','#FFCC0D');
c = hex2rgb(c);

%Proscillaridin A PK
nexttile
hold on
for i = 1:size(ProA.doses,2)
	plot(timeofDiv/60,ProA.TreatmentCc(i,:),'LineWidth',2,'Color',c(i,:))
end
hold off
xlim([0,72])
ylabel('Plasma Concentration (nM)')
set(gca,'yscale','log','FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('Proscillaridin A','FontSize',20)

%Digoxin PK
nexttile
hold on
for i = 1:size(Dig.doses,2)
    plot(timeofDiv/60,Dig.TreatmentCc(i,:),'LineWidth',2,'Color',c(i,:))
end
hold off
xlim([0,72])
set(gca,'yscale','log','FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('Digoxin','FontSize',20)

%Ouabain PK
nexttile
hold on
for i = 1:size(Oua.doses,2)
    plot(timeofDiv/60,Oua.TreatmentCc(i,:),'LineWidth',2,'Color',c(i,:))
end
hold off
legend(sprintf('%g nM',CarGly_doses_nM(1,1)),sprintf('%g nM',CarGly_doses_nM(1,2)),sprintf('%g nM',CarGly_doses_nM(1,3)),sprintf('%g nM',CarGly_doses_nM(1,4)),'Location', 'eastoutside','FontSize',17)
xlim([0,72])
set(gca,'yscale','log','FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('Ouabain','FontSize',20)

%Proscillaridin A toxicity
nexttile
hold on
for i = 1:size(ProA.doses,2)
    plot(timeofDiv/60,ProA.TreatmentToxicity(i,:),'LineWidth',2,'Color',c(i,:))
end
hold off
xlabel('Time (hours)')
xlim([0,72])
ylabel([{'Inhibition of K^+ Uptake'},{'in RBCs (%)'}])
ylim([0,110])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])

%Digoxin toxicity
nexttile
hold on
for i = 1:size(Dig.doses,2)
    plot(timeofDiv/60,Dig.TreatmentToxicity(i,:),'LineWidth',2,'Color',c(i,:))
end
hold off
xlabel('Time (hours)')
xlim([0,72])
ylim([0,110])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])

%Ouabain toxicity
nexttile
hold on
for i = 1:size(Oua.doses,2)
    plot(timeofDiv/60,Oua.TreatmentToxicity(i,:),'LineWidth',2,'Color',c(i,:))
end
hold off
xlabel('Time (hours)')
xlim([0,72])
ylim([0,110])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])

%% Figure 4B: Predicted PKPD responses of candidate glucocorticoids

figure
tiledlayout(1,2,'TileSpacing','compact');
c = char('#FF194D','#FFCC0D','#6CADA1','#2A5F65');
c = hex2rgb(c);

%Budesonide PK
nexttile
hold on
for i = 1:size(Bud.doses,2)
    plot(timeofDiv/60,Bud.TreatmentCc(i,:),'LineWidth',2,'Color',c(i,:))
end
hold off
xlabel('Time (hours)')
xlim([0,72])
ylabel('Plasma Concentration (nM)')
set(gca,'yscale','log','FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('Budesonide','FontSize',20)

%Mometasone PK
nexttile
hold on
for i = 1:size(Mom.doses,2)
    plot(timeofDiv/60,Mom.TreatmentCc(i,:),'LineWidth',2,'Color',c(i,:))
end
hold off
legend(sprintf('%g nM',Glu_doses_nM(1,1)),sprintf('%g nM',Glu_doses_nM(1,2)),sprintf('%g nM',Glu_doses_nM(1,3)),sprintf('%g nM',Glu_doses_nM(1,4)),'Location', 'eastoutside','FontSize',16)
xlabel('Time (hours)')
xlim([0,72])
set(gca,'yscale','log','FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('Mometasone','FontSize',20)
