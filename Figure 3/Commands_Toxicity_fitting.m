%% Fits of toxicity model for each cardiac glycoside

%load data
load('ProA_toxicity_data.mat');
load('Digoxin_toxicity_data.mat');
load('Ouabain_toxicity_data.mat');

ProA.units = (1/530.64968)*10^9; %mg/mL to nM
Dig.units = (1/780.938)*10^9;    %mg/mL to nM
Oua.units = (1/584.6525)*10^9;   %mg/mL to nM

modelfn_Tox = @(beta,C) beta(1)*(C.^beta(3))./(beta(2).^beta(3) + C.^beta(3));
nboot = 500; % bootstrap times

% set up data
ProA.Tox_dataConc = ProA_toxicity_data(:,1)*1e09; %in nM
ProA.Tox_data = ProA_toxicity_data(:,2);
ProA.Tox_dataUB = ProA_toxicity_data(:,3);
ProA.Tox_dataLB = ProA_toxicity_data(:,4);

Dig.Tox_dataConc = Digoxin_toxicity_data(:,1)*1e09; %in nM
Dig.Tox_data = Digoxin_toxicity_data(:,2);
Dig.Tox_dataUB = Digoxin_toxicity_data(:,3);
Dig.Tox_dataLB = Digoxin_toxicity_data(:,4);

Oua.Tox_dataConc = Ouabain_toxicity_data(:,1)*1e09; %in nM;
Oua.Tox_data = Ouabain_toxicity_data(:,2);
Oua.Tox_dataUB = Ouabain_toxicity_data(:,3);
Oua.Tox_dataLB = Ouabain_toxicity_data(:,4);

%% Parameter estimation

%Proscillaridin A
ProA.Tox_Emax_guess = 100;
ProA.Tox_IC50_guess = 9e-11;
ProA.Tox_h_guess = 2;

ProA.Tox_Para_guess = [ProA.Tox_Emax_guess ProA.Tox_IC50_guess ProA.Tox_h_guess]; 
[ProA.Tox_Para_fit,ProA.Tox_residual,ProA.Tox_J] = fitting_toxicity(ProA.Tox_dataConc,ProA.Tox_data,ProA.Tox_Para_guess);

ProA.Tox_Emax = ProA.Tox_Para_fit(1);
ProA.Tox_IC50 = ProA.Tox_Para_fit(2);
ProA.Tox_h    = ProA.Tox_Para_fit(3);

%Digoxin
Dig.Tox_Emax_guess = 100;
Dig.Tox_IC50_guess = 2e-9; 
Dig.Tox_h_guess = 0.89;

Dig.Tox_Para_guess = [Dig.Tox_Emax_guess Dig.Tox_IC50_guess Dig.Tox_h_guess]; 
[Dig.Tox_Para_fit,Dig.Tox_residual,Dig.Tox_J] = fitting_toxicity(Dig.Tox_dataConc,Dig.Tox_data,Dig.Tox_Para_guess);

Dig.Tox_Emax = Dig.Tox_Para_fit(1);
Dig.Tox_IC50 = Dig.Tox_Para_fit(2);
Dig.Tox_h    = Dig.Tox_Para_fit(3);

%Ouabain
Oua.Tox_Emax_guess = 100;
Oua.Tox_IC50_guess = 6e-10;
Oua.Tox_h_guess = 1.5;

Oua.Tox_Para_guess = [Oua.Tox_Emax_guess Oua.Tox_IC50_guess Oua.Tox_h_guess]; 
[Oua.Tox_Para_fit,Oua.Tox_residual,Oua.Tox_J] = fitting_toxicity(Oua.Tox_dataConc,Oua.Tox_data,Oua.Tox_Para_guess);

Oua.Tox_Emax = Oua.Tox_Para_fit(1);
Oua.Tox_IC50 = Oua.Tox_Para_fit(2);
Oua.Tox_h    = Oua.Tox_Para_fit(3);

%% Simulating toxicity model

%Proscillaridin A
ProA.Tox_simConc = logspace(-3,5,5000);
ProA.Tox_fit = modelfn_Tox(ProA.Tox_Para_fit,ProA.Tox_simConc);
ProA.Tox_fitdata = modelfn_Tox(ProA.Tox_Para_fit,ProA.Tox_dataConc);

%Digoxin
Dig.Tox_simConc = logspace(-3,5,1000);
Dig.Tox_fit = modelfn_Tox(Dig.Tox_Para_fit,Dig.Tox_simConc);
Dig.Tox_fitdata = modelfn_Tox(Dig.Tox_Para_fit,Dig.Tox_dataConc);

%Ouabain
Oua.Tox_simConc = logspace(-3,5,1000);
Oua.Tox_fit = modelfn_Tox(Oua.Tox_Para_fit,Oua.Tox_simConc);
Oua.Tox_fitdata = modelfn_Tox(Oua.Tox_Para_fit,Oua.Tox_dataConc);

%% Bootstrap of residuals: 95% confidence bands of fit

%Proscillaridin A
bootfn_Tox = @(conc,response)fitting_toxicity(conc,response,ProA.Tox_Para_guess);
[~,ProA.Tox_indicesb] = bootstrp(nboot,bootfn_Tox,ProA.Tox_dataConc',ProA.Tox_data');

ProA.Tox_residsb = ProA.Tox_residual(ProA.Tox_indicesb');
ProA.Toxb = repmat(ProA.Tox_fitdata',nboot,1) + ProA.Tox_residsb;

ProA.Tox_Parab = zeros(nboot,3);
ProA.Toxb_sim = zeros(nboot,size(ProA.Tox_simConc,2));
for j = 1:nboot
   [ProA.Tox_Parab(j,:),~,~] = fitting_toxicity(ProA.Tox_dataConc',ProA.Toxb(j,:),ProA.Tox_Para_guess);
   ProA.Toxb_sim(j,:) = modelfn_Tox(ProA.Tox_Parab(j,:),ProA.Tox_simConc);
end

ProA.Tox_cib = prctile(ProA.Tox_Parab,[2.5 97.5]); % Parameters' 95% bootstrap confidence intervals

ProA.Toxbci_sim = zeros(size(ProA.Tox_simConc,2),2); % For 95% bootstrap confidence bands
for i = 1:size(ProA.Tox_simConc,2)
   ProA.Toxbci_sim(i,:) = prctile(ProA.Toxb_sim(:,i),[2.5,97.5]); 
end

%Digoxin
bootfn_Tox = @(conc,response)fitting_toxicity(conc,response,Dig.Tox_Para_guess);
[~,Dig.Tox_indicesb] = bootstrp(nboot,bootfn_Tox,Dig.Tox_dataConc',Dig.Tox_data');

Dig.Tox_residsb = Dig.Tox_residual(Dig.Tox_indicesb');
Dig.Toxb = repmat(Dig.Tox_fitdata',nboot,1) + Dig.Tox_residsb;

Dig.Tox_Parab = zeros(nboot,3);
Dig.Toxb_sim = zeros(nboot,size(Dig.Tox_simConc,2));
for j = 1:nboot
   [Dig.Tox_Parab(j,:),~,~] = fitting_toxicity(Dig.Tox_dataConc',Dig.Toxb(j,:),Dig.Tox_Para_guess);
   Dig.Toxb_sim(j,:) = modelfn_Tox(Dig.Tox_Parab(j,:),Dig.Tox_simConc);
end

Dig.Tox_cib = prctile(Dig.Tox_Parab,[2.5 97.5]); % Parameters' 95% bootstrap confidence intervals

Dig.Toxbci_sim = zeros(size(Dig.Tox_simConc,2),2); % For 95% bootstrap confidence bands
for i = 1:size(Dig.Tox_simConc,2)
   Dig.Toxbci_sim(i,:) = prctile(Dig.Toxb_sim(:,i),[2.5,97.5]); 
end

%Ouabain
bootfn_Tox = @(conc,response)fitting_toxicity(conc,response,Oua.Tox_Para_guess);
[~,Oua.Tox_indicesb] = bootstrp(nboot,bootfn_Tox,Oua.Tox_dataConc',Oua.Tox_data');

Oua.Tox_residsb = Oua.Tox_residual(Oua.Tox_indicesb');
Oua.Toxb = repmat(Oua.Tox_fitdata',nboot,1) + Oua.Tox_residsb;

Oua.Tox_Parab = zeros(nboot,3);
Oua.Toxb_sim = zeros(nboot,size(Oua.Tox_simConc,2));
for j = 1:nboot
   [Oua.Tox_Parab(j,:),~,~] = fitting_toxicity(Oua.Tox_dataConc',Oua.Toxb(j,:),Oua.Tox_Para_guess);
   Oua.Toxb_sim(j,:) = modelfn_Tox(Oua.Tox_Parab(j,:),Oua.Tox_simConc);
end

Oua.Tox_cib = prctile(Oua.Tox_Parab,[2.5 97.5]); % Parameters' 95% bootstrap confidence intervals

Oua.Toxbci_sim = zeros(size(Oua.Tox_simConc,2),2); % For 95% bootstrap confidence bands
for i = 1:size(Oua.Tox_simConc,2)
   Oua.Toxbci_sim(i,:) = prctile(Oua.Toxb_sim(:,i),[2.5,97.5]); 
end

%% Save work

% save('ProA_toxicity.mat', '-struct', 'ProA');
% save('Digoxin_toxicity.mat', '-struct', 'Dig');
% save('Ouabain_toxicity.mat', '-struct', 'Oua');

%% Figure 3: Dose-dependent toxicity of cardiac glycosides

figure
tiledlayout(1,3,'TileSpacing','compact');

%Proscillaridin A
nexttile
hold on 
errorbar(ProA.Tox_dataConc,ProA.Tox_data,ProA.Tox_dataLB,ProA.Tox_dataUB,'o','MarkerEdgeColor','#001253','MarkerFaceColor','#001253','MarkerSize',8,'LineWidth',1.0,'Color','#001253'); %plotting the data;
plot(ProA.Tox_simConc,ProA.Tox_fit,'Color','#001253','LineWidth',1.5); % plotting curve
patch([ProA.Tox_simConc, fliplr(ProA.Tox_simConc)],[ProA.Toxbci_sim(:,1)', fliplr(ProA.Toxbci_sim(:,2)')],1,'facecolor','#001253', 'edgecolor', 'none', 'facealpha', 0.1); %CI
hold off
set(gca,'xscale','log','FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
set(gca, 'ytick', 0:20:100);
xlim([1e-3 1e+5])
ylim([-5 110])
xlabel('Concentration (nM)')
ylabel('Toxicity (% Inhibition of RBC K^+)')
title('Proscillaridin A','FontSize',20)

%Digoxin
nexttile
hold on 
errorbar(Dig.Tox_dataConc,Dig.Tox_data,Dig.Tox_dataLB,Dig.Tox_dataUB,'o','MarkerEdgeColor','#001253','MarkerFaceColor','#001253','MarkerSize',8,'LineWidth',1.0,'Color','#001253'); %plotting the data;
plot(Dig.Tox_simConc,Dig.Tox_fit,'Color','#001253','LineWidth',1.5); % plotting curve
patch([Dig.Tox_simConc,fliplr(Dig.Tox_simConc)],[Dig.Toxbci_sim(:,1)',fliplr(Dig.Toxbci_sim(:,2)')],1,'facecolor','#001253', 'edgecolor', 'none', 'facealpha', 0.1); %CI
hold off
set(gca,'xscale','log','FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
set(gca, 'ytick', 0:20:100);
xlim([1e-3 1e+5])
ylim([-5 110])
xlabel('Concentration (nM)')
ylabel('Toxicity (% Inhibition of RBC K^+)')
title('Digoxin','FontSize',20)

%Ouabain
nexttile
hold on 
errorbar(Oua.Tox_dataConc,Oua.Tox_data,Oua.Tox_dataLB,Oua.Tox_dataUB,'o','MarkerEdgeColor','#001253','MarkerFaceColor','#001253','MarkerSize',8,'LineWidth',1.0,'Color','#001253'); 
plot(Oua.Tox_simConc,Oua.Tox_fit,'Color','#001253','LineWidth',1.5); % plotting curve
patch([Oua.Tox_simConc,fliplr(Oua.Tox_simConc)],[Oua.Toxbci_sim(:,1)',fliplr(Oua.Toxbci_sim(:,2)')],1,'facecolor','#001253', 'edgecolor', 'none', 'facealpha', 0.1); %CI
hold off
set(gca,'xscale','log','FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
set(gca, 'ytick', 0:20:100);
xlim([1e-3 1e+5])
ylim([-5 110])
xlabel('Concentration (nM)')
ylabel('Toxicity (% Inhibition of RBC K^+)')
title('Ouabain','FontSize',20)