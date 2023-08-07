%% Fits of cell viability model for each anti-LSC drug

%load data
load('ProA_cell_viability_data.mat');
load('Digoxin_cell_viability_data.mat');
load('Ouabain_cell_viability_data.mat');
load('Budesonide_cell_viability_data.mat');
load('Mometasone_cell_viability_data.mat');

ProA.units = (1/530.64968)*10^9; %mg/mL to nM
Dig.units = (1/780.938)*10^9;    %mg/mL to nM
Oua.units = (1/584.6525)*10^9;   %mg/mL to nM
Bud.units = (1/430.534)*10^9;    %mg/mL to nM
Mom.units = (1/521.4)*10^9;      %mg/mL to nM

modelfn_Viab = @(beta,C) beta(2) - (beta(2)*C.^beta(4))./(beta(3).^beta(4) + C.^beta(4)) + beta(1);
nboot = 500; % bootstrap times

%set up data
ProA.v_dataConc = ProA_cell_viability_data(:,1)'; %in nM
ProA.vHSC_data = ProA_cell_viability_data(:,2)';
ProA.vHSC_dataUB = ProA_cell_viability_data(:,3)';
ProA.vHSC_dataLB = ProA_cell_viability_data(:,4)';
ProA.vLSC_data = ProA_cell_viability_data(:,5)'; 
ProA.vLSC_dataUB = ProA_cell_viability_data(:,6)';
ProA.vLSC_dataLB = ProA_cell_viability_data(:,7)';

Dig.v_dataConc = Digoxin_cell_viability_data(:,1)'; %in nM
Dig.vHSC_data = Digoxin_cell_viability_data(:,2)';
Dig.vHSC_dataUB = Digoxin_cell_viability_data(:,3)';
Dig.vHSC_dataLB = Digoxin_cell_viability_data(:,4)';
Dig.vLSC_data = Digoxin_cell_viability_data(:,5)'; 
Dig.vLSC_dataUB = Digoxin_cell_viability_data(:,6)';
Dig.vLSC_dataLB = Digoxin_cell_viability_data(:,7)';

Oua.v_dataConc = Ouabain_cell_viability_data(:,1)'; %in nM
Oua.vHSC_data = Ouabain_cell_viability_data(:,2)';
Oua.vHSC_dataUB = Ouabain_cell_viability_data(:,3)';
Oua.vHSC_dataLB = Ouabain_cell_viability_data(:,4)';
Oua.vLSC_data = Ouabain_cell_viability_data(:,5)'; 
Oua.vLSC_dataUB = Ouabain_cell_viability_data(:,6)';
Oua.vLSC_dataLB = Ouabain_cell_viability_data(:,7)';

Bud.v_dataConc = Budesonide_cell_viability_data(:,1)'; %in nM  
Bud.vHSC_data = Budesonide_cell_viability_data(:,2)';
Bud.vHSC_dataUB = Budesonide_cell_viability_data(:,3)';
Bud.vHSC_dataLB = Budesonide_cell_viability_data(:,4)';
Bud.vLSC_data = Budesonide_cell_viability_data(:,5)'; 
Bud.vLSC_dataUB = Budesonide_cell_viability_data(:,6)';
Bud.vLSC_dataLB = Budesonide_cell_viability_data(:,7)';

Mom.v_dataConc = Mometasone_cell_viability_data(:,1)'; %in nM 
Mom.vHSC_data = Mometasone_cell_viability_data(:,2)'; 
Mom.vHSC_dataUB = Mometasone_cell_viability_data(:,3)';
Mom.vHSC_dataLB = Mometasone_cell_viability_data(:,4)';
Mom.vLSC_data = Mometasone_cell_viability_data(:,5)'; 
Mom.vLSC_dataUB = Mometasone_cell_viability_data(:,6)';
Mom.vLSC_dataLB = Mometasone_cell_viability_data(:,7)';

%% Parameter estimation

%Proscillaridin A
ProA.vHSC_Emin_guess = 0;
ProA.vHSC_Emax_guess = ProA_cell_viability_data(1,2);
ProA.vHSC_IC50_guess = 29;
ProA.vHSC_h_guess = 2;
ProA.vLSC_Emin_guess = 0;
ProA.vLSC_Emax_guess = ProA_cell_viability_data(1,5);
ProA.vLSC_IC50_guess = 15;
ProA.vLSC_h_guess = 2;

ProA.vHSC_Para_guess = [ProA.vHSC_Emin_guess ProA.vHSC_Emax_guess ProA.vHSC_IC50_guess ProA.vHSC_h_guess]; 
ProA.vLSC_Para_guess = [ProA.vLSC_Emin_guess ProA.vLSC_Emax_guess ProA.vLSC_IC50_guess ProA.vLSC_h_guess]; 
[ProA.vHSC_Para_fit,ProA.vHSC_residual,ProA.vHSC_J] = fitting_IC50_curve(ProA.v_dataConc,ProA.vHSC_data,ProA.vHSC_Para_guess);
[ProA.vLSC_Para_fit,ProA.vLSC_residual,ProA.vLSC_J] = fitting_IC50_curve(ProA.v_dataConc,ProA.vLSC_data,ProA.vLSC_Para_guess);

ProA.vHSC_Emin = ProA.vHSC_Para_fit(1);
ProA.vHSC_Emax = ProA.vHSC_Para_fit(2);
ProA.vHSC_IC50 = ProA.vHSC_Para_fit(3);
ProA.vHSC_h    = ProA.vHSC_Para_fit(4);
ProA.vLSC_Emin = ProA.vLSC_Para_fit(1);
ProA.vLSC_Emax = ProA.vLSC_Para_fit(2);
ProA.vLSC_IC50 = ProA.vLSC_Para_fit(3);
ProA.vLSC_h    = ProA.vLSC_Para_fit(4);

ProA.vHSC_p = size(ProA.vHSC_Para_fit,2); % Number of parameters
ProA.vLSC_p = size(ProA.vLSC_Para_fit,2);

%Digoxin
Dig.vHSC_Emin_guess = 0;
Dig.vHSC_Emax_guess = Digoxin_cell_viability_data(1,2);
Dig.vHSC_IC50_guess = 43.40;
Dig.vHSC_h_guess = 2;
Dig.vLSC_Emin_guess = 0;
Dig.vLSC_Emax_guess = Digoxin_cell_viability_data(1,5);
Dig.vLSC_IC50_guess = 29.57;
Dig.vLSC_h_guess = 2;

Dig.vHSC_Para_guess = [Dig.vHSC_Emin_guess Dig.vHSC_Emax_guess Dig.vHSC_IC50_guess Dig.vHSC_h_guess];
Dig.vLSC_Para_guess = [Dig.vLSC_Emin_guess Dig.vLSC_Emax_guess Dig.vLSC_IC50_guess Dig.vLSC_h_guess];
[Dig.vHSC_Para_fit,Dig.vHSC_residual,Dig.vHSC_J] = fitting_IC50_curve(Dig.v_dataConc,Dig.vHSC_data,Dig.vHSC_Para_guess);
[Dig.vLSC_Para_fit,Dig.vLSC_residual,Dig.vLSC_J] = fitting_IC50_curve(Dig.v_dataConc,Dig.vLSC_data,Dig.vLSC_Para_guess);

Dig.vHSC_Emin = Dig.vHSC_Para_fit(1);
Dig.vHSC_Emax = Dig.vHSC_Para_fit(2);
Dig.vHSC_IC50 = Dig.vHSC_Para_fit(3);
Dig.vHSC_h    = Dig.vHSC_Para_fit(4);
Dig.vLSC_Emin = Dig.vLSC_Para_fit(1);
Dig.vLSC_Emax = Dig.vLSC_Para_fit(2);
Dig.vLSC_IC50 = Dig.vLSC_Para_fit(3);
Dig.vLSC_h    = Dig.vLSC_Para_fit(4);

Dig.vHSC_p = size(Dig.vHSC_Para_fit,2); % Number of parameters
Dig.vLSC_p = size(Dig.vLSC_Para_fit,2);

%Ouabain
Oua.vHSC_Emin_guess = 0;
Oua.vHSC_Emax_guess = Ouabain_cell_viability_data(1,2);
Oua.vHSC_IC50_guess = 31.45;
Oua.vHSC_h_guess = 2;
Oua.vLSC_Emin_guess = 0;
Oua.vLSC_Emax_guess = Ouabain_cell_viability_data(1,5);
Oua.vLSC_IC50_guess = 15.17;
Oua.vLSC_h_guess = 2;

Oua.vHSC_Para_guess = [Oua.vHSC_Emin_guess Oua.vHSC_Emax_guess Oua.vHSC_IC50_guess Oua.vHSC_h_guess];
Oua.vLSC_Para_guess = [Oua.vLSC_Emin_guess Oua.vLSC_Emax_guess Oua.vLSC_IC50_guess Oua.vLSC_h_guess];
[Oua.vHSC_Para_fit,Oua.vHSC_residual,Oua.vHSC_J] = fitting_IC50_curve(Oua.v_dataConc,Oua.vHSC_data,Oua.vHSC_Para_guess);
[Oua.vLSC_Para_fit,Oua.vLSC_residual,Oua.vLSC_J] = fitting_IC50_curve(Oua.v_dataConc,Oua.vLSC_data,Oua.vLSC_Para_guess);

Oua.vHSC_Emin = Oua.vHSC_Para_fit(1);
Oua.vHSC_Emax = Oua.vHSC_Para_fit(2);
Oua.vHSC_IC50 = Oua.vHSC_Para_fit(3);
Oua.vHSC_h    = Oua.vHSC_Para_fit(4);
Oua.vLSC_Emin = Oua.vLSC_Para_fit(1);
Oua.vLSC_Emax = Oua.vLSC_Para_fit(2);
Oua.vLSC_IC50 = Oua.vLSC_Para_fit(3);
Oua.vLSC_h    = Oua.vLSC_Para_fit(4);

Oua.vHSC_p = size(Oua.vHSC_Para_fit,2); % Number of parameters
Oua.vLSC_p = size(Oua.vLSC_Para_fit,2);

%Budesonide
Bud.vHSC_Emin_guess = 0;
Bud.vHSC_Emax_guess = Budesonide_cell_viability_data(1,2);
Bud.vHSC_IC50_guess = 69.95;
Bud.vHSC_h_guess = 2;
Bud.vLSC_Emin_guess = 0;
Bud.vLSC_Emax_guess = Budesonide_cell_viability_data(1,5);
Bud.vLSC_IC50_guess = 4.04;
Bud.vLSC_h_guess = 2;

Bud.vHSC_Para_guess = [Bud.vHSC_Emin_guess Bud.vHSC_Emax_guess Bud.vHSC_IC50_guess Bud.vHSC_h_guess];
Bud.vLSC_Para_guess = [Bud.vLSC_Emin_guess Bud.vLSC_Emax_guess Bud.vLSC_IC50_guess Bud.vLSC_h_guess];
[Bud.vHSC_Para_fit,Bud.vHSC_residual,Bud.vHSC_J] = fitting_IC50_curve(Bud.v_dataConc,Bud.vHSC_data,Bud.vHSC_Para_guess);
[Bud.vLSC_Para_fit,Bud.vLSC_residual,Bud.vLSC_J] = fitting_IC50_curve(Bud.v_dataConc,Bud.vLSC_data,Bud.vLSC_Para_guess);

Bud.vHSC_Emin = Bud.vHSC_Para_fit(1);
Bud.vHSC_Emax = Bud.vHSC_Para_fit(2);
Bud.vHSC_IC50 = Bud.vHSC_Para_fit(3);
Bud.vHSC_h    = Bud.vHSC_Para_fit(4);
Bud.vLSC_Emin = Bud.vLSC_Para_fit(1);
Bud.vLSC_Emax = Bud.vLSC_Para_fit(2);
Bud.vLSC_IC50 = Bud.vLSC_Para_fit(3);
Bud.vLSC_h    = Bud.vLSC_Para_fit(4);

Bud.vHSC_p = size(Bud.vHSC_Para_fit,2); % Number of parameters
Bud.vLSC_p = size(Bud.vLSC_Para_fit,2);

%Mometasone
Mom.vHSC_Emin_guess = 0;
Mom.vHSC_Emax_guess = Mometasone_cell_viability_data(1,2);
Mom.vHSC_IC50_guess = 50.76;
Mom.vHSC_h_guess = 2;
Mom.vLSC_Emin_guess = 0;
Mom.vLSC_Emax_guess = Mometasone_cell_viability_data(1,5);
Mom.vLSC_IC50_guess = 0.86;
Mom.vLSC_h_guess = 2;

Mom.vHSC_Para_guess = [Mom.vHSC_Emin_guess Mom.vHSC_Emax_guess Mom.vHSC_IC50_guess Mom.vHSC_h_guess];
Mom.vLSC_Para_guess = [Mom.vLSC_Emin_guess Mom.vLSC_Emax_guess Mom.vLSC_IC50_guess Mom.vLSC_h_guess];
[Mom.vHSC_Para_fit,Mom.vHSC_residual,Mom.vHSC_J] = fitting_IC50_curve(Mom.v_dataConc,Mom.vHSC_data,Mom.vHSC_Para_guess);
[Mom.vLSC_Para_fit,Mom.vLSC_residual,Mom.vLSC_J] = fitting_IC50_curve(Mom.v_dataConc,Mom.vLSC_data,Mom.vLSC_Para_guess);

Mom.vHSC_Emin = Mom.vHSC_Para_fit(1);
Mom.vHSC_Emax = Mom.vHSC_Para_fit(2);
Mom.vHSC_IC50 = Mom.vHSC_Para_fit(3);
Mom.vHSC_h    = Mom.vHSC_Para_fit(4);
Mom.vLSC_Emin = Mom.vLSC_Para_fit(1);
Mom.vLSC_Emax = Mom.vLSC_Para_fit(2);
Mom.vLSC_IC50 = Mom.vLSC_Para_fit(3);
Mom.vLSC_h    = Mom.vLSC_Para_fit(4);

Mom.vHSC_p = size(Mom.vHSC_Para_fit,2); % Number of parameters
Mom.vLSC_p = size(Mom.vLSC_Para_fit,2);

%% Simulating cell viability model

%Proscillaridin A
ProA.vHSC_simConc = linspace(0,100,1000);  
ProA.vHSC_fit = modelfn_Viab(ProA.vHSC_Para_fit,ProA.vHSC_simConc); 
ProA.vHSC_datafit = modelfn_Viab(ProA.vHSC_Para_fit,ProA.v_dataConc);
ProA.vLSC_simConc = linspace(0,100,1000);
ProA.vLSC_fit = modelfn_Viab(ProA.vLSC_Para_fit,ProA.vLSC_simConc); 
ProA.vLSC_datafit = modelfn_Viab(ProA.vLSC_Para_fit,ProA.v_dataConc);

%Digoxin
Dig.vHSC_simConc = linspace(0,100,1000);
Dig.vHSC_fit = modelfn_Viab(Dig.vHSC_Para_fit,Dig.vHSC_simConc); 
Dig.vHSC_datafit = modelfn_Viab(Dig.vHSC_Para_fit,Dig.v_dataConc);
Dig.vLSC_simConc = linspace(0,100,1000);  
Dig.vLSC_fit = modelfn_Viab(Dig.vLSC_Para_fit,Dig.vLSC_simConc); 
Dig.vLSC_datafit = modelfn_Viab(Dig.vLSC_Para_fit,Dig.v_dataConc);

%Ouabain
Oua.vHSC_simConc = linspace(0,100,1000);
Oua.vHSC_fit = modelfn_Viab(Oua.vHSC_Para_fit,Oua.vHSC_simConc); 
Oua.vHSC_datafit = modelfn_Viab(Oua.vHSC_Para_fit,Oua.v_dataConc);
Oua.vLSC_simConc = linspace(0,100,1000);
Oua.vLSC_fit = modelfn_Viab(Oua.vLSC_Para_fit,Oua.vLSC_simConc);
Oua.vLSC_datafit = modelfn_Viab(Oua.vLSC_Para_fit,Oua.v_dataConc);

%Budesonide
Bud.vHSC_simConc = logspace(-2,3,1000); 
Bud.vHSC_fit = modelfn_Viab(Bud.vHSC_Para_fit,Bud.vHSC_simConc);
Bud.vHSC_datafit = modelfn_Viab(Bud.vHSC_Para_fit,Bud.v_dataConc);
Bud.vLSC_simConc = logspace(-2,3,1000); 
Bud.vLSC_fit = modelfn_Viab(Bud.vLSC_Para_fit,Bud.vLSC_simConc);
Bud.vLSC_datafit = modelfn_Viab(Bud.vLSC_Para_fit,Bud.v_dataConc);

%Mometasone
Mom.vHSC_simConc = logspace(-2,3,1000);
Mom.vHSC_fit = modelfn_Viab(Mom.vHSC_Para_fit,Mom.vHSC_simConc);
Mom.vHSC_datafit = modelfn_Viab(Mom.vHSC_Para_fit,Mom.v_dataConc);
Mom.vLSC_simConc = logspace(-2,3,1000);
Mom.vLSC_fit = modelfn_Viab(Mom.vLSC_Para_fit,Mom.vLSC_simConc);
Mom.vLSC_datafit = modelfn_Viab(Mom.vLSC_Para_fit,Mom.v_dataConc);

%% Bootstrap of residuals: 95% confidence bands of fit

%Proscillaridin A
bootfn_Viab = @(conc,response)fitting_IC50_curve(conc,response,ProA.vHSC_Para_guess);
[~,ProA.vHSC_indicesb] = bootstrp(nboot,bootfn_Viab,ProA.v_dataConc,ProA.vHSC_data); 
bootfn_Viab = @(conc,response)fitting_IC50_curve(conc,response,ProA.vLSC_Para_guess);
[~,ProA.vLSC_indicesb] = bootstrp(nboot,bootfn_Viab,ProA.v_dataConc,ProA.vLSC_data); 

ProA.vHSC_residsb = ProA.vHSC_residual(ProA.vHSC_indicesb');
ProA.vLSC_residsb = ProA.vLSC_residual(ProA.vLSC_indicesb');
ProA.vHSCb = repmat(ProA.vHSC_datafit,nboot,1) + ProA.vHSC_residsb;
ProA.vLSCb = repmat(ProA.vLSC_datafit,nboot,1) + ProA.vLSC_residsb;

ProA.vHSC_Parab = zeros(nboot,ProA.vHSC_p);
ProA.vLSC_Parab = zeros(nboot,ProA.vLSC_p);
ProA.vHSCb_sim = zeros(nboot,size(ProA.vHSC_simConc,2));
ProA.vLSCb_sim = zeros(nboot,size(ProA.vLSC_simConc,2));
for j = 1:nboot
   [ProA.vHSC_Parab(j,:),~,~] = fitting_IC50_curve(ProA.v_dataConc,ProA.vHSCb(j,:),ProA.vHSC_Para_guess);
   [ProA.vLSC_Parab(j,:),~,~] = fitting_IC50_curve(ProA.v_dataConc,ProA.vLSCb(j,:),ProA.vLSC_Para_guess);
   ProA.vHSCb_sim(j,:) = modelfn_Viab(ProA.vHSC_Parab(j,:),ProA.vHSC_simConc); 
   ProA.vLSCb_sim(j,:) = modelfn_Viab(ProA.vLSC_Parab(j,:),ProA.vLSC_simConc);
end

ProA.vHSC_cib = prctile(ProA.vHSC_Parab,[2.5 97.5]); % Parameters' 95% bootstrap confidence intervals
ProA.vLSC_cib = prctile(ProA.vLSC_Parab,[2.5 97.5]);

ProA.vHSCbci_sim = zeros(size(ProA.vHSC_simConc,2),2); % For 95% bootstrap confidence bands
ProA.vLSCbci_sim = zeros(size(ProA.vLSC_simConc,2),2);
for i = 1:size(ProA.vHSC_simConc,2)
   ProA.vHSCbci_sim(i,:) = prctile(ProA.vHSCb_sim(:,i),[2.5,97.5]);
   ProA.vLSCbci_sim(i,:) = prctile(ProA.vLSCb_sim(:,i),[2.5,97.5]);
end

%Digoxin
bootfn_Viab = @(conc,response)fitting_IC50_curve(conc,response,Dig.vHSC_Para_guess);
[~,Dig.vHSC_indicesb] = bootstrp(nboot,bootfn_Viab,Dig.v_dataConc,Dig.vHSC_data);
bootfn_Viab = @(conc,response)fitting_IC50_curve(conc,response,Dig.vLSC_Para_guess);
[~,Dig.vLSC_indicesb] = bootstrp(nboot,bootfn_Viab,Dig.v_dataConc,Dig.vLSC_data);

Dig.vHSC_residsb = Dig.vHSC_residual(Dig.vHSC_indicesb');
Dig.vLSC_residsb = Dig.vLSC_residual(Dig.vLSC_indicesb');
Dig.vHSCb = repmat(Dig.vHSC_datafit,nboot,1) + Dig.vHSC_residsb;
Dig.vLSCb = repmat(Dig.vLSC_datafit,nboot,1) + Dig.vLSC_residsb;

Dig.vHSC_Parab = zeros(nboot,Dig.vHSC_p);
Dig.vLSC_Parab = zeros(nboot,Dig.vLSC_p);
Dig.vHSCb_sim = zeros(nboot,size(Dig.vHSC_simConc,2));
Dig.vLSCb_sim = zeros(nboot,size(Dig.vLSC_simConc,2));
for j = 1:nboot
   [Dig.vHSC_Parab(j,:),~,~] = fitting_IC50_curve(Dig.v_dataConc,Dig.vHSCb(j,:),Dig.vHSC_Para_guess);
   [Dig.vLSC_Parab(j,:),~,~] = fitting_IC50_curve(Dig.v_dataConc,Dig.vLSCb(j,:),Dig.vLSC_Para_guess);
   Dig.vHSCb_sim(j,:) = modelfn_Viab(Dig.vHSC_Parab(j,:),Dig.vHSC_simConc); 
   Dig.vLSCb_sim(j,:) = modelfn_Viab(Dig.vLSC_Parab(j,:),Dig.vLSC_simConc);
end

Dig.vHSC_cib = prctile(Dig.vHSC_Parab,[2.5 97.5]); % Parameters' 95% bootstrap confidence intervals
Dig.vLSC_cib = prctile(Dig.vLSC_Parab,[2.5 97.5]);

Dig.vHSCbci_sim = zeros(size(Dig.vHSC_simConc,2),2); % For 95% bootstrap confidence bands
Dig.vLSCbci_sim = zeros(size(Dig.vLSC_simConc,2),2);
for i = 1:size(Dig.vHSC_simConc,2)
   Dig.vHSCbci_sim(i,:) = prctile(Dig.vHSCb_sim(:,i),[2.5,97.5]);
   Dig.vLSCbci_sim(i,:) = prctile(Dig.vLSCb_sim(:,i),[2.5,97.5]);
end

%Ouabain
bootfn_Viab = @(conc,response)fitting_IC50_curve(conc,response,Oua.vHSC_Para_guess);
[~,Oua.vHSC_indicesb] = bootstrp(nboot,bootfn_Viab,Oua.v_dataConc,Oua.vHSC_data);
bootfn_Viab = @(conc,response)fitting_IC50_curve(conc,response,Oua.vLSC_Para_guess);
[~,Oua.vLSC_indicesb] = bootstrp(nboot,bootfn_Viab,Oua.v_dataConc,Oua.vLSC_data);

Oua.vHSC_residsb = Oua.vHSC_residual(Oua.vHSC_indicesb');
Oua.vLSC_residsb = Oua.vLSC_residual(Oua.vLSC_indicesb');
Oua.vHSCb = repmat(Oua.vHSC_datafit,nboot,1) + Oua.vHSC_residsb;
Oua.vLSCb = repmat(Oua.vLSC_datafit,nboot,1) + Oua.vLSC_residsb;

Oua.vHSC_Parab = zeros(nboot,Oua.vHSC_p);
Oua.vLSC_Parab = zeros(nboot,Oua.vLSC_p);
Oua.vHSCb_sim = zeros(nboot,size(Oua.vHSC_simConc,2));
Oua.vLSCb_sim = zeros(nboot,size(Oua.vLSC_simConc,2));
for j = 1:nboot
   [Oua.vHSC_Parab(j,:),~,~] = fitting_IC50_curve(Oua.v_dataConc,Oua.vHSCb(j,:),Oua.vHSC_Para_guess);
   [Oua.vLSC_Parab(j,:),~,~] = fitting_IC50_curve(Oua.v_dataConc,Oua.vLSCb(j,:),Oua.vLSC_Para_guess);
   Oua.vHSCb_sim(j,:) = modelfn_Viab(Oua.vHSC_Parab(j,:),Oua.vHSC_simConc); 
   Oua.vLSCb_sim(j,:) = modelfn_Viab(Oua.vLSC_Parab(j,:),Oua.vLSC_simConc);
end

Oua.vHSC_cib = prctile(Oua.vHSC_Parab,[2.5 97.5]); % Parameters' 95% bootstrap confidence intervals
Oua.vLSC_cib = prctile(Oua.vLSC_Parab,[2.5 97.5]);

Oua.vHSCbci_sim = zeros(size(Oua.vHSC_simConc,2),2); % For 95% bootstrap confidence bands
Oua.vLSCbci_sim = zeros(size(Oua.vLSC_simConc,2),2);
for i = 1:size(Oua.vHSC_simConc,2)
   Oua.vHSCbci_sim(i,:) = prctile(Oua.vHSCb_sim(:,i),[2.5,97.5]);
   Oua.vLSCbci_sim(i,:) = prctile(Oua.vLSCb_sim(:,i),[2.5,97.5]);
end

%Budesonide
bootfn_Viab = @(conc,response)fitting_IC50_curve(conc,response,Bud.vHSC_Para_guess);
[~,Bud.vHSC_indicesb] = bootstrp(nboot,bootfn_Viab,Bud.v_dataConc,Bud.vHSC_data);
bootfn_Viab = @(conc,response)fitting_IC50_curve(conc,response,Bud.vLSC_Para_guess);
[~,Bud.vLSC_indicesb] = bootstrp(nboot,bootfn_Viab,Bud.v_dataConc,Bud.vLSC_data);

Bud.vHSC_residsb = Bud.vHSC_residual(Bud.vHSC_indicesb');
Bud.vLSC_residsb = Bud.vLSC_residual(Bud.vLSC_indicesb');
Bud.vHSCb = repmat(Bud.vHSC_datafit,nboot,1) + Bud.vHSC_residsb;
Bud.vLSCb = repmat(Bud.vLSC_datafit,nboot,1) + Bud.vLSC_residsb;

Bud.vHSC_Parab = zeros(nboot,Bud.vHSC_p);
Bud.vLSC_Parab = zeros(nboot,Bud.vLSC_p);
Bud.vHSCb_sim = zeros(nboot,size(Bud.vHSC_simConc,2));
Bud.vLSCb_sim = zeros(nboot,size(Bud.vLSC_simConc,2));
for j = 1:nboot
   [Bud.vHSC_Parab(j,:),~,~] = fitting_IC50_curve(Bud.v_dataConc,Bud.vHSCb(j,:),Bud.vHSC_Para_guess); 
   [Bud.vLSC_Parab(j,:),~,~] = fitting_IC50_curve(Bud.v_dataConc,Bud.vLSCb(j,:),Bud.vLSC_Para_guess);
   Bud.vHSCb_sim(j,:) = modelfn_Viab(Bud.vHSC_Parab(j,:),Bud.vHSC_simConc);
   Bud.vLSCb_sim(j,:) = modelfn_Viab(Bud.vLSC_Parab(j,:),Bud.vLSC_simConc);
end

Bud.vHSC_cib = prctile(Bud.vHSC_Parab,[2.5 97.5]); % Parameters' 95% bootstrap confidence intervals
Bud.vLSC_cib = prctile(Bud.vLSC_Parab,[2.5 97.5]);

Bud.vHSCbci_sim = zeros(size(Bud.vHSC_simConc,2),2); % For 95% bootstrap confidence bands
Bud.vLSCbci_sim = zeros(size(Bud.vLSC_simConc,2),2);
for i = 1:size(Bud.vHSC_simConc,2)
   Bud.vHSCbci_sim(i,:) = prctile(Bud.vHSCb_sim(:,i),[2.5,97.5]);
   Bud.vLSCbci_sim(i,:) = prctile(Bud.vLSCb_sim(:,i),[2.5,97.5]);
end

%Mometasone
bootfn_Viab = @(conc,response)fitting_IC50_curve(conc,response,Mom.vHSC_Para_guess);
[~,Mom.vHSC_indicesb] = bootstrp(nboot,bootfn_Viab,Mom.v_dataConc,Mom.vHSC_data);
bootfn_Viab = @(conc,response)fitting_IC50_curve(conc,response,Mom.vLSC_Para_guess);
[~,Mom.vLSC_indicesb] = bootstrp(nboot,bootfn_Viab,Mom.v_dataConc,Mom.vLSC_data);

Mom.vHSC_residsb = Mom.vHSC_residual(Mom.vHSC_indicesb');
Mom.vLSC_residsb = Mom.vLSC_residual(Mom.vLSC_indicesb');
Mom.vHSCb = repmat(Mom.vHSC_data,nboot,1) + Mom.vHSC_residsb;
Mom.vLSCb = repmat(Mom.vLSC_data,nboot,1) + Mom.vLSC_residsb;

Mom.vHSC_Parab = zeros(nboot,Mom.vHSC_p);
Mom.vLSC_Parab = zeros(nboot,Mom.vLSC_p);
Mom.vHSCb_sim = zeros(nboot,size(Mom.vHSC_simConc,2));
Mom.vLSCb_sim = zeros(nboot,size(Mom.vLSC_simConc,2));
for j = 1:nboot
   [Mom.vHSC_Parab(j,:),~,~] = fitting_IC50_curve(Mom.v_dataConc,Mom.vHSCb(j,:),Mom.vHSC_Para_guess);
   [Mom.vLSC_Parab(j,:),~,~] = fitting_IC50_curve(Mom.v_dataConc,Mom.vLSCb(j,:),Mom.vLSC_Para_guess);
   Mom.vHSCb_sim(j,:) = modelfn_Viab(Mom.vHSC_Parab(j,:),Mom.vHSC_simConc);
   Mom.vLSCb_sim(j,:) = modelfn_Viab(Mom.vLSC_Parab(j,:),Mom.vLSC_simConc);
end

Mom.vHSC_cib = prctile(Mom.vHSC_Parab,[2.5 97.5]); % Parameters' 95% bootstrap confidence intervals
Mom.vLSC_cib = prctile(Mom.vLSC_Parab,[2.5 97.5]);

Mom.vHSCbci_sim = zeros(size(Mom.vHSC_simConc,2),2); % For 95% bootstrap confidence bands
Mom.vLSCbci_sim = zeros(size(Mom.vLSC_simConc,2),2);
for i = 1:size(Mom.vHSC_simConc,2)
   Mom.vHSCbci_sim(i,:) = prctile(Mom.vHSCb_sim(:,i),[2.5,97.5]);
   Mom.vLSCbci_sim(i,:) = prctile(Mom.vLSCb_sim(:,i),[2.5,97.5]);
end

%% Simulating fitness advantage

%Proscillaridin A
ProA.sd = 1-ProA.vLSC_fit./ProA.vHSC_fit - (1-ProA.vLSC_fit(1,1)./ProA.vHSC_fit(1,1));
ProA.sd = max(ProA.sd,0);
ProA.sp = 1.25*(ProA.vLSC_fit./ProA.vLSC_fit(1,1));

%Digoxin
Dig.sd = 1-Dig.vLSC_fit./Dig.vHSC_fit - (1-Dig.vLSC_fit(1,1)./Dig.vHSC_fit(1,1));
Dig.sp = 1.25*(Dig.vLSC_fit./Dig.vLSC_fit(1,1));

%Ouabain
Oua.sd = 1-Oua.vLSC_fit./Oua.vHSC_fit  - (1-Oua.vLSC_fit(1,1)./Oua.vHSC_fit(1,1));
Oua.sp = 1.25*(Oua.vLSC_fit./Oua.vLSC_fit(1,1));

%Budesonide
Bud.sd = 1-Bud.vLSC_fit./Bud.vHSC_fit - (1-Bud.vLSC_fit(1,1)./Bud.vHSC_fit(1,1));
Bud.sp = 1.25*(Bud.vLSC_fit./Bud.vLSC_fit(1,1));

%Mometasone
Mom.sd = 1-Mom.vLSC_fit./Mom.vHSC_fit - (1-Mom.vLSC_fit(1,1)./Mom.vHSC_fit(1,1));
Mom.sp = 1.25*(Mom.vLSC_fit./Mom.vLSC_fit(1,1));

%% Save work

% save('ProA_cell_viability.mat', '-struct', 'ProA');
% save('Digoxin_cell_viability.mat', '-struct', 'Dig');
% save('Ouabain_cell_viability.mat', '-struct', 'Oua');
% save('Budesonide_cell_viability.mat', '-struct', 'Bud');
% save('Mometasone_cell_viability.mat', '-struct', 'Mom');

%% Figure 2A: Cell Viability

figure
tiledlayout(2,3,'TileSpacing','compact');

%Proscillaridin A
nexttile 
hold on 
h(1) = errorbar(ProA.v_dataConc,ProA.vHSC_data,ProA.vHSC_dataLB,ProA.vHSC_dataUB,'^','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',1.0,'Color','b'); %plotting the data
h(2) = plot(ProA.vHSC_simConc,ProA.vHSC_fit,'b','LineWidth',1.5); % plotting curve
h(3) = patch([ProA.vHSC_simConc, fliplr(ProA.vHSC_simConc)],[ProA.vHSCbci_sim(:,1)', fliplr(ProA.vHSCbci_sim(:,2)')],1,'facecolor', 'b', 'edgecolor', 'none', 'facealpha', 0.1); %CI
h(4) = errorbar(ProA.v_dataConc,ProA.vLSC_data,ProA.vLSC_dataLB,ProA.vLSC_dataUB,'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10,'LineWidth',1.0,'Color','r'); %plotting the data
h(5) = plot(ProA.vLSC_simConc,ProA.vLSC_fit,'r','LineWidth',1.5); % plotting curve
h(6) = patch([ProA.vLSC_simConc, fliplr(ProA.vLSC_simConc)],[ProA.vLSCbci_sim(:,1)', fliplr(ProA.vLSCbci_sim(:,2)')],1,'facecolor', 'r', 'edgecolor', 'none', 'facealpha', 0.1); %CI
hold off
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
xlim([0 100])
ylim([-5 120])
xlabel('Concentration (nM)')
ylabel('Cell Viability (%)')
title('Proscillaridin A','FontSize',20)

%Digoxin
nexttile
hold on  
h(1) = errorbar(Dig.v_dataConc,Dig.vHSC_data,Dig.vHSC_dataLB,Dig.vHSC_dataUB,'^','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',1.0,'Color','b'); %plotting the data
h(2) = plot(Dig.vHSC_simConc,Dig.vHSC_fit,'b','LineWidth',1.5); % plotting curve
h(3) = patch([Dig.vHSC_simConc, fliplr(Dig.vHSC_simConc)],[Dig.vHSCbci_sim(:,1)', fliplr(Dig.vHSCbci_sim(:,2)')],1,'facecolor', 'b', 'edgecolor', 'none', 'facealpha', 0.1); %CI
h(4) = errorbar(Dig.v_dataConc,Dig.vLSC_data,Dig.vLSC_dataLB,Dig.vLSC_dataUB,'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10,'LineWidth',1.0,'Color','r'); %plotting the data
h(5) = plot(Dig.vLSC_simConc,Dig.vLSC_fit,'r','LineWidth',1.5); % plotting curve
h(6) = patch([Dig.vLSC_simConc, fliplr(Dig.vLSC_simConc)],[Dig.vLSCbci_sim(:,1)', fliplr(Dig.vLSCbci_sim(:,2)')],1,'facecolor', 'r', 'edgecolor', 'none', 'facealpha', 0.1); %CI
hold off
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
xlim([0 100])
ylim([-5 120])
xlabel('Concentration (nM)')
ylabel('Cell Viability (%)')
title('Digoxin','FontSize',20)

%Ouabain
nexttile 
hold on  
h(1) = errorbar(Oua.v_dataConc,Oua.vHSC_data,Oua.vHSC_dataLB,Oua.vHSC_dataUB,'^','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',1.0,'Color','b'); %plotting the data
h(2) = plot(Oua.vHSC_simConc,Oua.vHSC_fit,'b','LineWidth',1.5); % plotting curve
h(3) = patch([Oua.vHSC_simConc, fliplr(Oua.vHSC_simConc)],[Oua.vHSCbci_sim(:,1)', fliplr(Oua.vHSCbci_sim(:,2)')],1,'facecolor', 'b', 'edgecolor', 'none', 'facealpha', 0.1); %CI
h(4) = errorbar(Oua.v_dataConc,Oua.vLSC_data,Oua.vLSC_dataLB,Oua.vLSC_dataUB,'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10,'LineWidth',1.0,'Color','r'); %plotting the data
h(5) = plot(Oua.vLSC_simConc,Oua.vLSC_fit,'r','LineWidth',1.5); % plotting curve
h(6) = patch([Oua.vLSC_simConc, fliplr(Oua.vLSC_simConc)],[Oua.vLSCbci_sim(:,1)', fliplr(Oua.vLSCbci_sim(:,2)')],1,'facecolor', 'r', 'edgecolor', 'none', 'facealpha', 0.1); %CI
hold off
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
set(gca, 'ytick', 0:20:130);
xlim([0 100])
ylim([-5 130])
xlabel('Concentration (nM)')
ylabel('Cell Viability (%)')
title('Ouabain','FontSize',20)

%Budesonide
nexttile 
hold on  
h(1) = errorbar(Bud.v_dataConc,Bud.vHSC_data,Bud.vHSC_dataLB,Bud.vHSC_dataUB,'^','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',1.0,'Color','b'); %plotting the data
h(2) = plot(Bud.vHSC_simConc,Bud.vHSC_fit,'b','LineWidth',1.5); % plotting curve
h(3) = patch([Bud.vHSC_simConc,fliplr(Bud.vHSC_simConc)],[Bud.vHSCbci_sim(:,1)',fliplr(Bud.vHSCbci_sim(:,2)')],1,'facecolor', 'b', 'edgecolor', 'none', 'facealpha', 0.1); %CI
h(4) = errorbar(Bud.v_dataConc,Bud.vLSC_data,Bud.vLSC_dataLB,Bud.vLSC_dataUB,'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10,'LineWidth',1.0,'Color','r'); %plotting the data
h(5) = plot(Bud.vLSC_simConc,Bud.vLSC_fit,'r','LineWidth',1.5); % plotting curve
h(6) = patch([Bud.vLSC_simConc,fliplr(Bud.vLSC_simConc)],[Bud.vLSCbci_sim(:,1)',fliplr(Bud.vLSCbci_sim(:,2)')],1,'facecolor', 'r', 'edgecolor', 'none', 'facealpha', 0.1); %CI
hold off
set(gca,'xscale','log','FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
xlim([0 1000])
ylim([-5 120])
xlabel('Concentration (nM)')
ylabel('Cell Viability (%)')
title('Budesonide','FontSize',20)

%Mometasone
nexttile 
hold on  
h(1) = errorbar(Mom.v_dataConc,Mom.vHSC_data,Mom.vHSC_dataLB,Mom.vHSC_dataUB,'^','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',1.0,'Color','b'); %plotting the data
h(2) = plot(Mom.vHSC_simConc,Mom.vHSC_fit,'b','LineWidth',1.5); % plotting curve
h(3) = patch([Mom.vHSC_simConc,fliplr(Mom.vHSC_simConc)],[Mom.vHSCbci_sim(:,1)',fliplr(Mom.vHSCbci_sim(:,2)')],1,'facecolor', 'b', 'edgecolor', 'none', 'facealpha', 0.1); %CI
h(4) = errorbar(Mom.v_dataConc,Mom.vLSC_data,Mom.vLSC_dataLB,Mom.vLSC_dataUB,'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10,'LineWidth',1.0,'Color','r'); %plotting the data
h(5) = plot(Mom.vLSC_simConc,Mom.vLSC_fit,'r','LineWidth',1.5); % plotting curve
h(6) = patch([Mom.vLSC_simConc,fliplr(Mom.vLSC_simConc)],[Mom.vLSCbci_sim(:,1)',fliplr(Mom.vLSCbci_sim(:,2)')],1,'facecolor', 'r', 'edgecolor', 'none', 'facealpha', 0.1); %CI

l(1) = plot(nan,nan,'^','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',10); % plotting curve
l(2) = plot(nan,nan,'s','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',10); % plotting curve
l(3) = plot(nan,nan,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',10); % plotting curve
l(4) = plot(nan,nan,'b','LineWidth',1.5); % plotting curve
l(5) = plot(nan,nan,'r','LineWidth',1.5); % plotting curve
hold off
legend(l, {'CD34- AML 8227','CD34+ Cord Blood','CD34+ AML 8227','Fit for HSC', 'Fit for LSC'},'Orientation','horizontal','Location', 'southoutside','FontSize',14)
set(gca,'xscale','log','FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
xlim([0 1000])
ylim([-5 120])
xlabel('Concentration (nM)')
ylabel('Cell Viability (%)')
title('Mometasone','FontSize',20)

%% Figure 2B: Fitness Advantages

c = char('#0B84A5','#F6C85F');
c = hex2rgb(c);

figure
tiledlayout(2,3,'TileSpacing','compact');

%Proscillaridin A
nexttile 
hold on
plot(ProA.vHSC_simConc,ProA.sd,'LineWidth',2,'Color',c(2,:))
plot(ProA.vHSC_simConc,ProA.sp,'LineWidth',2,'Color',c(1,:))
hold off
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
set(gca, 'ytick', 0:0.25:1.5);
xlim([0 100])
ylim([-0.1 1.50])
xlabel('Concentration (nM)')
ylabel('Fitness Value')
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('Proscillaridin A','FontSize',20)

%Digoxin
nexttile 
hold on
plot(Dig.vHSC_simConc,Dig.sd,'LineWidth',2,'Color',c(2,:))
plot(Dig.vHSC_simConc,Dig.sp,'LineWidth',2,'Color',c(1,:))
hold off
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
set(gca, 'ytick', 0:0.25:1.5);
xlim([0 100])
ylim([-0.1 1.50])
xlabel('Concentration (nM)')
ylabel('Fitness Value')
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('Digoxin','FontSize',20)

%Ouabain
nexttile 
hold on
plot(Oua.vHSC_simConc,Oua.sd,'LineWidth',2,'Color',c(2,:));
plot(Oua.vHSC_simConc,Oua.sp,'LineWidth',2,'Color',c(1,:));
hold off
set(gca, 'ytick', 0:0.25:1.5);
xlim([0 100])
ylim([-0.1 1.50])
xlabel('Concentration (nM)')
ylabel('Fitness Value')
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('Ouabain','FontSize',20)

%Budesonide
nexttile 
hold on
plot(Bud.vHSC_simConc,Bud.sd,'LineWidth',2,'Color',c(2,:))
plot(Bud.vHSC_simConc,Bud.sp,'LineWidth',2,'Color',c(1,:))
hold off
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
set(gca, 'ytick', 0:0.25:1.5);
xlim([0 1000])
ylim([-0.1 1.50])
xlabel('Concentration (nM)')
ylabel('Fitness Value')
set(gca,'xscale','log','FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('Budesonide','FontSize',20)

%Mometasone
nexttile 
hold on
legend_plot1 = plot(Mom.vHSC_simConc,Mom.sd,'LineWidth',2,'Color',c(2,:));
legend_plot2 = plot(Mom.vHSC_simConc,Mom.sp,'LineWidth',2,'Color',c(1,:));
hold off
legend([legend_plot1 legend_plot2], {'s_d','s_p'},'Orientation','vertical','Location','eastoutside','FontSize',14)
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
set(gca, 'ytick', 0:0.25:1.5);
xlim([0 1000])
ylim([-0.1 1.50])
xlabel('Concentration (nM)')
ylabel('Fitness Value')
set(gca,'xscale','log','FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('Mometasone','FontSize',20)
