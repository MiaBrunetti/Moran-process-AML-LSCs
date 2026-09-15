%% Commands for Figures

%load work form Commands_Treatment_Model.mat
ProA = load('ProA_Moran.mat');
Dig = load('Dig_Moran.mat');
Oua = load('Oua_Moran.mat');
Bud = load('Bud_Moran.mat');
Mom = load('Mom_Moran.mat');
AraC = load('AraC_Moran.mat');

%parameters of Moran process
sim_num = 100;                                %number of simulation of the Moran Process
sim_time = 365*24*60;                         %in min
Nd = 50000;                                   %total number of stem cells
r = 28*24*60;                                 %division rate of one stem cell in minutes: 1 HCS divides every r minutes (28 days)
div = round((1/(r/Nd))*100)/100;              %number of cell divisions in a minute
frac = 0.08;                                  %initial fraction of LSCs
maxNoEvents = sim_time*div;                   %number of events/divisions in the Moran Process for the simulated time
timeofDiv = linspace(0,sim_time,maxNoEvents); %time at which each division occurs for a full year

%setting doses 
AraC_doses_nM = [125 25 5 1]; %in nM
CarGly_doses_nM = [50 30 20 10]; %in nM
Glu_doses_nM = [25 10 1.5 0.25]; %in nM

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

%% Figure 2A: Cell Viability

F2A = figure;
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
ylim([-5 130])
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
ylim([-5 130])
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
ylim([-5 130])
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
hold off
set(gca,'xscale','log','FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
xlim([0 1000])
ylim([-5 130])
xlabel('Concentration (nM)')
ylabel('Cell Viability (%)')
title('Mometasone','FontSize',20)

nexttile
hold on 
%h(1) = scatter(AraC.v_dataConc,AraC.vHSC_data,'o','MarkerEdgeColor','k','MarkerFaceColor','#8EB1DC');
h(1) = errorbar(AraC.v_meanConc,AraC.vHSC_mean,AraC.vHSC_dataLB,AraC.vHSC_dataUB,'^','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',1.0,'Color','b'); %plotting the data
h(2) = plot(AraC.vHSC_simConc,AraC.vHSC_fit,'b','LineWidth',1.5); % plotting curve
h(3) = patch([AraC.vHSC_simConc,fliplr(AraC.vHSC_simConc)],[AraC.vHSCbci_sim(:,1)',fliplr(AraC.vHSCbci_sim(:,2)')],1,'facecolor', 'b', 'edgecolor', 'none', 'facealpha', 0.1); %CI
%h(5) = scatter(AraC.v_dataConc,AraC.vLSC_data,'o','MarkerEdgeColor','k','MarkerFaceColor','#E58A8C');
h(4) = errorbar(AraC.v_meanConc,AraC.vLSC_mean,AraC.vLSC_dataLB,AraC.vLSC_dataUB,'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10,'LineWidth',1.0,'Color','r'); %plotting the data
h(5) = plot(AraC.vLSC_simConc,AraC.vLSC_fit,'r','LineWidth',1.5); % plotting curve
h(6) = patch([AraC.vLSC_simConc,fliplr(AraC.vLSC_simConc)],[AraC.vLSCbci_sim(:,1)',fliplr(AraC.vLSCbci_sim(:,2)')],1,'facecolor', 'r', 'edgecolor', 'none', 'facealpha', 0.1); %CI

l(1) = plot(nan,nan,'^','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',10); % plotting curve
l(2) = plot(nan,nan,'s','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',10); % plotting curve
l(3) = plot(nan,nan,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',10); % plotting curve
l(4) = plot(nan,nan,'b','LineWidth',1.5); % plotting curve
l(5) = plot(nan,nan,'r','LineWidth',1.5); % plotting curve
hold off
% legend(l, {'CD34- AML 8227','CD34+ Cord Blood','CD34+ AML 8227','Fit for HSC', 'Fit for LSC'},'Orientation','horizontal','Location', 'southoutside','FontSize',14)
set(gca,'xscale','log','FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
xlim([1e-02 1000])
ylim([-5 130])
xlabel('Concentration (nM)')
ylabel('Cell Viability (%)')
title('Cytarabine','FontSize',20)

%% Figure 2B: Fitness Advantages

c = char('#0B84A5','#F6C85F');
c = hex2rgb(c);

F2B = figure;
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
plot(Mom.vHSC_simConc,Mom.sd,'LineWidth',2,'Color',c(2,:));
plot(Mom.vHSC_simConc,Mom.sp,'LineWidth',2,'Color',c(1,:));
hold off
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
set(gca, 'ytick', 0:0.25:1.5);
xlim([0 1000])
ylim([-0.1 1.50])
xlabel('Concentration (nM)')
ylabel('Fitness Value')
set(gca,'xscale','log','FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('Mometasone','FontSize',20)

%Cytarabine
nexttile
hold on
legend_plot1 = plot(AraC.vHSC_simConc,max(0,AraC.sd),'LineWidth',2,'Color',c(2,:));
legend_plot2 = plot(AraC.vHSC_simConc,AraC.sp,'LineWidth',2,'Color',c(1,:));
hold off
% legend([legend_plot1 legend_plot2], {'s_d','s_p'},'Orientation','vertical','Location','eastoutside','FontSize',14)
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
set(gca, 'ytick', 0:0.25:1.5);
xlim([0 1000])
ylim([-0.1 1.50])
xlabel('Concentration (nM)')
ylabel('Fitness Value')
set(gca,'xscale','log','FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('Cytarabine','FontSize',20)

%% Figure 4: Predicted PKPD responses

F4A = figure;
tiledlayout(2,3,'TileSpacing','compact');
c_CarGly = char('#702A8C','#BF2669','#FF7326','#FFCC0D');
c_CarGly = hex2rgb(c_CarGly);

%Proscillaridin A PK
nexttile
hold on
for i = 1:size(ProA.doses,2)
	plot(timeofDiv/60,ProA.TreatmentCc(i,:),'LineWidth',2,'Color',c_CarGly(i,:))
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
    plot(timeofDiv/60,Dig.TreatmentCc(i,:),'LineWidth',2,'Color',c_CarGly(i,:))
end
hold off
xlim([0,72])
set(gca,'yscale','log','FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('Digoxin','FontSize',20)

%Ouabain PK
nexttile
hold on
for i = 1:size(Oua.doses,2)
    plot(timeofDiv/60,Oua.TreatmentCc(i,:),'LineWidth',2,'Color',c_CarGly(i,:))
end
hold off
%legend(sprintf('%g nM',CarGly_doses_nM(1,1)),sprintf('%g nM',CarGly_doses_nM(1,2)),sprintf('%g nM',CarGly_doses_nM(1,3)),sprintf('%g nM',CarGly_doses_nM(1,4)),'Location', 'eastoutside','FontSize',17)
xlim([0,72])
set(gca,'yscale','log','FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('Ouabain','FontSize',20)

%Proscillaridin A toxicity
nexttile
hold on
for i = 1:size(ProA.doses,2)
    plot(timeofDiv/60,ProA.TreatmentToxicity(i,:),'LineWidth',2,'Color',c_CarGly(i,:))
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
    plot(timeofDiv/60,Dig.TreatmentToxicity(i,:),'LineWidth',2,'Color',c_CarGly(i,:))
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
    plot(timeofDiv/60,Oua.TreatmentToxicity(i,:),'LineWidth',2,'Color',c_CarGly(i,:))
end
hold off
xlabel('Time (hours)')
xlim([0,72])
ylim([0,110])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])

F4B = figure;
tiledlayout(1,3,'TileSpacing','compact');
c_Glu = char('#FF194D','#FFCC0D','#6CADA1','#2A5F65');
c_Glu = hex2rgb(c_Glu);

%Budesonide PK
nexttile
hold on
for i = 1:size(Bud.doses,2)
    plot(timeofDiv/60,Bud.TreatmentCc(i,:),'LineWidth',2,'Color',c_Glu(i,:))
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
    plot(timeofDiv/60,Mom.TreatmentCc(i,:),'LineWidth',2,'Color',c_Glu(i,:))
end
hold off
%legend(sprintf('%g nM',Glu_doses_nM(1,1)),sprintf('%g nM',Glu_doses_nM(1,2)),sprintf('%g nM',Glu_doses_nM(1,3)),sprintf('%g nM',Glu_doses_nM(1,4)),'Location', 'eastoutside','FontSize',16)
xlabel('Time (hours)')
xlim([0,72])
set(gca,'yscale','log','FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('Mometasone','FontSize',20)

%Cytarabine PK
c_AraC = char('#B33DC6','#27AEEF','#87BC45','#EF9B20');
c_AraC = hex2rgb(c_AraC);

nexttile
hold on
for i = 1:size(AraC.doses,2)
    plot(timeofDiv/60,AraC.TreatmentCc(i,:),'LineWidth',2,'Color',c_AraC(i,:))
end
hold off
%legend(sprintf('%g nM',AraC_doses_nM(1,1)),sprintf('%g nM',AraC_doses_nM(1,2)),sprintf('%g nM',AraC_doses_nM(1,3)),sprintf('%g nM',AraC_doses_nM(1,4)),'Location', 'eastoutside','FontSize',16)
xlabel('Time (hours)')
ylabel('Plasma Concentration (nM)')
xlim([0,72])
ylim([1e-05,2e+02])
set(gca,'yscale','log','FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('Cytarabine','FontSize',20)

%% Figure 5A: LSC expansion and viability of cardiac glycodides

c_CarGly = char('#702A8C','#BF2669','#FF7326','#FFCC0D','#000000');
color_CarGly = hex2rgb(c_CarGly);
options.color_area = color_CarGly;
options.color_line = color_CarGly;
options.legend = [CarGly_doses_nM 0];
options.alpha = 0.1;
options.line_width = 2;
options.error = 'std';
options.x_axis = timeofDiv/(60*24);

F5A = figure;
tiledlayout(2,3,'TileSpacing','compact');

nexttile %Proscillaridin A
plot_areaerrorbar_multiple(ProA.LSC,options)
ylabel('LSC Number')
xlim([0 sim_time/(60*24)])
ylim([0 35000])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('Proscillaridin A','FontSize',20)

nexttile %Digoxin
plot_areaerrorbar_multiple(Dig.LSC,options)
xlim([0 sim_time/(60*24)])
ylim([0 35000])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('Digoxin','FontSize',20)

nexttile %Ouabain
plot_areaerrorbar_multiple(Oua.LSC,options)
legend('Orientation','vertical','Location', 'eastoutside','FontSize',18)
xlim([0 sim_time/(60*24)])
ylim([0 35000])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('Ouabain','FontSize',20)

nexttile %Proscillaridin A
plot_areaerrorbar_multiple(ProA.vLSC,options)
xlabel('Time (days)')
ylabel('LSC Viability (%)')
xlim([0 sim_time/(60*24)])
ylim([0 120])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])

nexttile %Digoxin
plot_areaerrorbar_multiple(Dig.vLSC,options)
xlabel('Time (days)')
xlim([0 sim_time/(60*24)])
ylim([0 120])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])

nexttile %Ouabain
plot_areaerrorbar_multiple(Oua.vLSC,options)
xlabel('Time (days)')
xlim([0 sim_time/(60*24)])
ylim([0 120])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])

%% Figure 5B-C: LSC expansion and viability of glucocorticoids and cytarabine

c_Glu = char('#FF194D','#FFCC0D','#6CADA1','#2A5F65','#000000');
color_Glu = hex2rgb(c_Glu);
options.legend = [Glu_doses_nM 0];
options.color_area = color_Glu;
options.color_line = color_Glu;
options.alpha = 0.1;
options.line_width = 2;
options.error = 'std';
options.x_axis = timeofDiv/(60*24);

F5B = figure;
tiledlayout(2,3,'TileSpacing','compact');

nexttile %Budesonide
plot_areaerrorbar_multiple(Bud.LSC,options)
ylabel('LSC Number')
xlim([0 sim_time/(60*24)])
ylim([0 35000])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('Budesonide','FontSize',20)

nexttile %Mometasone
plot_areaerrorbar_multiple(Mom.LSC,options)
legend('Orientation','vertical','Location', 'eastoutside','FontSize',18)
xlim([0 sim_time/(60*24)])
ylim([0 35000])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('Mometasone','FontSize',20)

c_AraC = char('#B33DC6','#27AEEF','#87BC45','#EF9B20','#000000');
color_AraC = hex2rgb(c_AraC);
options.color_area = color_AraC;
options.color_line = color_AraC;
options.legend = [AraC_doses_nM 0];

nexttile %AraC
plot_areaerrorbar_multiple(AraC.LSC,options)
legend('Orientation','vertical','Location', 'eastoutside','FontSize',18)
ylabel('LSC Number')
xlim([0 sim_time/(60*24)])
ylim([0 35000])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('Cytarabine','FontSize',20)

options.color_area = color_Glu;
options.color_line = color_Glu;

nexttile %Budesonide
plot_areaerrorbar_multiple(Bud.vLSC,options)
xlabel('Time (days)')
ylabel('LSC Viability (%)')
xlim([0 sim_time/(60*24)])
ylim([0 120])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])

nexttile %Mometasone
plot_areaerrorbar_multiple(Mom.vLSC,options)
xlabel('Time (days)')
xlim([0 sim_time/(60*24)])
ylim([0 120])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])

options.color_area = color_AraC;
options.color_line = color_AraC;

nexttile %AraC
plot_areaerrorbar_multiple(AraC.vLSC,options)
xlabel('Time (days)')
ylabel('LSC Viability (%)')
xlim([0 sim_time/(60*24)])
ylim([0 120])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])

%% Figure 5D-G: Comparaison glucocorticoids to cytarabine response

% Mean value and standard deviation of Moran process at day 90, 180 & 365
index365 = 365*24*60*div;
index180 = 180*24*60*div;
index90 = 90*24*60*div;

AraC_Moran = AraC.LSC;
Bud_Moran = Bud.LSC;
Mom_Moran = Mom.LSC;

%Compare HSC# of lowest efficient doses of budesonide (10nM) and mometasone (1.5nM) vs Ara-C (25 and 125nM)
AraC_HSC365 = zeros(100,size(AraC_Moran,1));
AraC_HSC180 = zeros(100,size(AraC_Moran,1));
AraC_HSC90 = zeros(100,size(AraC_Moran,1));

Bud_HSC365 = zeros(100,size(Bud_Moran,1));
Bud_HSC180 = zeros(100,size(Bud_Moran,1));
Bud_HSC90 = zeros(100,size(Bud_Moran,1));

Mom_HSC365 = zeros(100,size(Mom_Moran,1));
Mom_HSC180 = zeros(100,size(Mom_Moran,1));
Mom_HSC90 = zeros(100,size(Mom_Moran,1));

for i = 1:size(AraC_Moran,1) %last index is no treatment
    data_LSC = AraC_Moran{i,1};

    % Finding LSC distributions at day 365, 183, and 90
    AraC_HSC365(:,i) = Nd - data_LSC(:,index365);
    AraC_HSC180(:,i) = Nd - data_LSC(:,index180);
    AraC_HSC90(:,i) = Nd - data_LSC(:,index90);
end

for i = 1:size(Bud_Moran,1) %last index is no treatment
    data_LSC = Bud_Moran{i,1};

    % Finding LSC distributions at day 365, 183, and 90
    Bud_HSC365(:,i) = Nd - data_LSC(:,index365);
    Bud_HSC180(:,i) = Nd - data_LSC(:,index180);
    Bud_HSC90(:,i) = Nd - data_LSC(:,index90);
end

for i = 1:size(Mom_Moran,1) %last index is no treatment
    data_LSC = Mom_Moran{i,1};

    % Finding LSC distributions at day 365, 183, and 90
    Mom_HSC365(:,i) = Nd - data_LSC(:,index365);
    Mom_HSC180(:,i) = Nd - data_LSC(:,index180);
    Mom_HSC90(:,i) = Nd - data_LSC(:,index90);
end

HSC_AraC125nM_AraC125nM = [AraC_HSC90(:,1)/mean(AraC_HSC90(:,1)) AraC_HSC180(:,1)/mean(AraC_HSC180(:,1)) AraC_HSC365(:,1)/mean(AraC_HSC365(:,1))];
HSC_AraC125nM_Bud10nM = [Bud_HSC90(:,2)/mean(AraC_HSC90(:,1)) Bud_HSC180(:,2)/mean(AraC_HSC180(:,1)) Bud_HSC365(:,2)/mean(AraC_HSC365(:,1))];
HSC_AraC125nM_Mom10nM = [Mom_HSC90(:,2)/mean(AraC_HSC90(:,1)) Mom_HSC180(:,2)/mean(AraC_HSC180(:,1)) Mom_HSC365(:,2)/mean(AraC_HSC365(:,1))];
HSCAraC125nM_Mom1_25nM = [Mom_HSC90(:,3)/mean(AraC_HSC90(:,1)) Mom_HSC180(:,3)/mean(AraC_HSC180(:,1)) Mom_HSC365(:,3)/mean(AraC_HSC365(:,1))];

HSC_AraC25nM_AraC25nM = [AraC_HSC90(:,2)/mean(AraC_HSC90(:,2)) AraC_HSC180(:,2)/mean(AraC_HSC180(:,2)) AraC_HSC365(:,2)/mean(AraC_HSC365(:,2))];
HSC_AraC25nM_Bud10nM = [Bud_HSC90(:,2)/mean(AraC_HSC90(:,2)) Bud_HSC180(:,2)/mean(AraC_HSC180(:,2)) Bud_HSC365(:,2)/mean(AraC_HSC365(:,2))];
HSC_AraC25nM_Mom10nM = [Mom_HSC90(:,2)/mean(AraC_HSC90(:,2)) Mom_HSC180(:,2)/mean(AraC_HSC180(:,2)) Mom_HSC365(:,2)/mean(AraC_HSC365(:,2))];
HSC_AraC25nM_Mom1_25nM = [Mom_HSC90(:,3)/mean(AraC_HSC90(:,2)) Mom_HSC180(:,3)/mean(AraC_HSC180(:,2)) Mom_HSC365(:,3)/mean(AraC_HSC365(:,2))];

HSC_y125 = [mean(HSC_AraC125nM_AraC125nM); mean(HSC_AraC125nM_Bud10nM); mean(HSCAraC125nM_Mom1_25nM)]';
HSC_err125 = [std(HSC_AraC125nM_AraC125nM); std(HSC_AraC125nM_Bud10nM); std(HSCAraC125nM_Mom1_25nM)]';
HSC_y25 = [mean(HSC_AraC25nM_AraC25nM); mean(HSC_AraC25nM_Bud10nM); mean(HSC_AraC25nM_Mom1_25nM)]';
HSC_err25 = [std(HSC_AraC25nM_AraC25nM); std(HSC_AraC25nM_Bud10nM); std(HSC_AraC25nM_Mom1_25nM)]';

bh_HSC125 = figure; 
bar(HSC_y125(:,3),'FaceColor',hex2rgb('#33a02c'),'LineWidth',1.5)
hold on 
bar(HSC_y125(:,2),'FaceColor',hex2rgb('#b2df8a'),'LineWidth',1.5)
bar(HSC_y125(:,1),'FaceColor',hex2rgb('#1f78b4'),'LineWidth',1.5)
errorbar(HSC_y125,HSC_err125,'.k','LineWidth',1.5);
bh(1) = bar(nan,nan,'FaceColor',hex2rgb('#a6cee3'));
bh(2) = bar(nan,nan,'FaceColor',hex2rgb('#1f78b4'));
bh(3) = bar(nan,nan,'FaceColor',hex2rgb('#b2df8a'));
bh(4) = bar(nan,nan,'FaceColor',hex2rgb('#33a02c'));
hold off
xticklabels({'90','180','365'})
legend(bh,({'25 nM Ara-C','125 nM Ara-C','10 nM Budesonide','1.5 nM Mometasone'}),'Location','bestoutside');
xlabel('Time (days)')
ylabel('Fold Change in HSC Number')
set(gca,'FontSize',18,'TickLength',[0.02 0.025])

bh_HSC25 = figure;
bar(HSC_y25(:,3),'FaceColor',hex2rgb('#33a02c'),'LineWidth',1.5)
hold on 
bar(HSC_y25(:,2),'FaceColor',hex2rgb('#b2df8a'),'LineWidth',1.5)
bar(HSC_y25(:,1),'FaceColor',hex2rgb('#a6cee3'),'LineWidth',1.5)
errorbar(HSC_y25,HSC_err25,'.k','LineWidth',1.5);
bh(1) = bar(nan,nan,'FaceColor',hex2rgb('#a6cee3'));
bh(2) = bar(nan,nan,'FaceColor',hex2rgb('#1f78b4'));
bh(3) = bar(nan,nan,'FaceColor',hex2rgb('#b2df8a'));
bh(4) = bar(nan,nan,'FaceColor',hex2rgb('#33a02c'));
hold off
xticklabels({'90','180','365'})
legend(bh,({'25 nM Ara-C','125 nM Ara-C','10 nM Budesonide','1.5 nM Mometasone'}),'Location','bestoutside');
xlabel('Time (days)')
ylabel('Fold Change in HSC Number')
set(gca,'FontSize',18,'TickLength',[0.02 0.025])

%Compare LSC# of lowest efficient doses of Budesonide (10nM) and mometasone (1.5nM) vs Ara-C (25 and 125nM)
AraC_LSC365 = zeros(100,size(AraC_Moran,1));
AraC_LSC180 = zeros(100,size(AraC_Moran,1));
AraC_LSC90 = zeros(100,size(AraC_Moran,1));

Bud_LSC365 = zeros(100,size(Bud_Moran,1));
Bud_LSC180 = zeros(100,size(Bud_Moran,1));
Bud_LSC90 = zeros(100,size(Bud_Moran,1));

Mom_LSC365 = zeros(100,size(Mom_Moran,1));
Mom_LSC180 = zeros(100,size(Mom_Moran,1));
Mom_LSC90 = zeros(100,size(Mom_Moran,1));

for i = 1:size(AraC_Moran,1) %last index is no treatment
    data_LSC = AraC_Moran{i,1};

    % Finding LSC distributions at day 365, 183, and 90
    AraC_LSC365(:,i) = data_LSC(:,index365);
    AraC_LSC180(:,i) = data_LSC(:,index180);
    AraC_LSC90(:,i) = data_LSC(:,index90);
end

for i = 1:size(Bud_Moran,1) %last index is no treatment
    data_LSC = Bud_Moran{i,1};

    % Finding LSC distributions at day 365, 183, and 90
    Bud_LSC365(:,i) = data_LSC(:,index365);
    Bud_LSC180(:,i) = data_LSC(:,index180);
    Bud_LSC90(:,i) = data_LSC(:,index90);
end

for i = 1:size(Mom_Moran,1) %last index is no treatment
    data_LSC = Mom_Moran{i,1};

    % Finding LSC distributions at day 365, 183, and 90
    Mom_LSC365(:,i) = data_LSC(:,index365);
    Mom_LSC180(:,i) = data_LSC(:,index180);
    Mom_LSC90(:,i) = data_LSC(:,index90);
end

LSC_AraC125nM_AraC125nM = [AraC_LSC90(:,1)/mean(AraC_LSC90(:,1)) AraC_LSC180(:,1)/mean(AraC_LSC180(:,1)) AraC_LSC365(:,1)/mean(AraC_LSC365(:,1))];
LSC_AraC125nM_Bud10nM = [Bud_LSC90(:,2)/mean(AraC_LSC90(:,1)) Bud_LSC180(:,2)/mean(AraC_LSC180(:,1)) Bud_LSC365(:,2)/mean(AraC_LSC365(:,1))];
LSC_AraC125nM_Mom10nM = [Mom_LSC90(:,2)/mean(AraC_LSC90(:,1)) Mom_LSC180(:,2)/mean(AraC_LSC180(:,1)) Mom_LSC365(:,2)/mean(AraC_LSC365(:,1))];
LSCAraC125nM_Mom1_25nM = [Mom_LSC90(:,3)/mean(AraC_LSC90(:,1)) Mom_LSC180(:,3)/mean(AraC_LSC180(:,1)) Mom_LSC365(:,3)/mean(AraC_LSC365(:,1))];

LSC_AraC25nM_AraC25nM = [AraC_LSC90(:,2)/mean(AraC_LSC90(:,2)) AraC_LSC180(:,2)/mean(AraC_LSC180(:,2)) AraC_LSC365(:,2)/mean(AraC_LSC365(:,2))];
LSC_AraC25nM_Bud10nM = [Bud_LSC90(:,2)/mean(AraC_LSC90(:,2)) Bud_LSC180(:,2)/mean(AraC_LSC180(:,2)) Bud_LSC365(:,2)/mean(AraC_LSC365(:,2))];
LSC_AraC25nM_Mom10nM = [Mom_LSC90(:,2)/mean(AraC_LSC90(:,2)) Mom_LSC180(:,2)/mean(AraC_LSC180(:,2)) Mom_LSC365(:,2)/mean(AraC_LSC365(:,2))];
LSC_AraC25nM_Mom1_25nM = [Mom_LSC90(:,3)/mean(AraC_LSC90(:,2)) Mom_LSC180(:,3)/mean(AraC_LSC180(:,2)) Mom_LSC365(:,3)/mean(AraC_LSC365(:,2))];

LSC_y125 = [mean(LSC_AraC125nM_AraC125nM); mean(LSC_AraC125nM_Bud10nM); mean(LSCAraC125nM_Mom1_25nM)]';
LSC_err125 = [std(LSC_AraC125nM_AraC125nM); std(LSC_AraC125nM_Bud10nM); std(LSCAraC125nM_Mom1_25nM)]';
LSC_y25 = [mean(LSC_AraC25nM_AraC25nM); mean(LSC_AraC25nM_Bud10nM); mean(LSC_AraC25nM_Mom1_25nM)]';
LSC_err25 = [std(LSC_AraC25nM_AraC25nM); std(LSC_AraC25nM_Bud10nM); std(LSC_AraC25nM_Mom1_25nM)]';

bh_LSC125 = figure; 
bar(LSC_y125(:,1),'FaceColor',hex2rgb('#1f78b4'),'LineWidth',1.5)
hold on 
bar(LSC_y125(:,2),'FaceColor',hex2rgb('#b2df8a'),'LineWidth',1.5)
bar(LSC_y125(:,3),'FaceColor',hex2rgb('#33a02c'),'LineWidth',1.5)
errorbar(LSC_y125,LSC_err125,'.k','LineWidth',1.5);
bh(1) = bar(nan,nan,'FaceColor',hex2rgb('#a6cee3'));
bh(2) = bar(nan,nan,'FaceColor',hex2rgb('#1f78b4'));
bh(3) = bar(nan,nan,'FaceColor',hex2rgb('#b2df8a'));
bh(4) = bar(nan,nan,'FaceColor',hex2rgb('#33a02c'));
hold off
xticklabels({'90','180','365'})
legend(bh,({'25 nM Ara-C','125 nM Ara-C','10 nM Budesonide','1.5 nM Mometasone'}),'Location','bestoutside');
xlabel('Time (days)')
ylabel('Fold Change in LSC Number')
set(gca,'FontSize',18,'TickLength',[0.02 0.025])

bh_LSC25 = figure;
bar(LSC_y25(:,1),'FaceColor',hex2rgb('#a6cee3'),'LineWidth',1.5)
hold on 
bar(LSC_y25(:,2),'FaceColor',hex2rgb('#b2df8a'),'LineWidth',1.5)
bar(LSC_y25(:,3),'FaceColor',hex2rgb('#33a02c'),'LineWidth',1.5)
errorbar(LSC_y25,LSC_err25,'.k','LineWidth',1.5);
bh(1) = bar(nan,nan,'FaceColor',hex2rgb('#a6cee3'));
bh(2) = bar(nan,nan,'FaceColor',hex2rgb('#1f78b4'));
bh(3) = bar(nan,nan,'FaceColor',hex2rgb('#b2df8a'));
bh(4) = bar(nan,nan,'FaceColor',hex2rgb('#33a02c'));
hold off
xticklabels({'90','180','365'})
legend(bh,({'25 nM Ara-C','125 nM Ara-C','10 nM Budesonide','1.5 nM Mometasone'}),'Location','bestoutside');
xlabel('Time (days)')
ylabel('Fold Change in LSC Number')
set(gca,'FontSize',18,'TickLength',[0.02 0.025])

%% T-test compare budesonide and mometasone to Ara-C

%LSC
[h90_LSC_AraC125_Bud10,p90_LSC_AraC125_Bud10] = ttest2(AraC_LSC90(:,1),Bud_LSC90(:,2),'Vartype','unequal','Alpha',0.01);
[h180_LSC_AraC125_Bud10,p180_LSC_AraC125_Bud10] = ttest2(AraC_LSC180(:,1),Bud_LSC180(:,2),'Vartype','unequal','Alpha',0.01);
[h365_LSC_AraC125_Bud10,p365_LSC_AraC125_Bud10] = ttest2(AraC_LSC365(:,1),Bud_LSC365(:,2),'Vartype','unequal','Alpha',0.01);

[h90_LSC_AraC125_Mom1_25,p90_LSC_AraC125_Mom1_25] = ttest2(AraC_LSC90(:,1),Mom_LSC90(:,3),'Vartype','unequal','Alpha',0.01);
[h180_LSC_AraC125_Mom1_25,p180_LSC_AraC125_Mom1_25] = ttest2(AraC_LSC180(:,1),Mom_LSC180(:,3),'Vartype','unequal','Alpha',0.01);
[h365_LSC_AraC125_Mom1_25,p365_LSC_AraC125_Mom1_25] = ttest2(AraC_LSC365(:,1),Mom_LSC365(:,3),'Vartype','unequal','Alpha',0.01);

[h90_LSC_AraC25_Bud10,p90_LSC_AraC25_Bud10] = ttest2(AraC_LSC90(:,2),Bud_LSC90(:,2),'Vartype','unequal','Alpha',0.01);
[h180_LSC_AraC25_Bud10,p180_LSC_AraC25_Bud10] = ttest2(AraC_LSC180(:,2),Bud_LSC180(:,2),'Vartype','unequal','Alpha',0.01);
[h365_LSC_AraC25_Bud10,p365_LSC_AraC25_Bud10] = ttest2(AraC_LSC365(:,2),Bud_LSC365(:,2),'Vartype','unequal','Alpha',0.01);

[h90_LSC_AraC25_Mom1_25,p90_LSC_AraC25_Mom1_25] = ttest2(AraC_LSC90(:,2),Mom_LSC90(:,3),'Vartype','unequal','Alpha',0.01);
[h180_LSC_AraC25_Mom1_25,p180_LSC_AraC25_Mom1_25] = ttest2(AraC_LSC180(:,2),Mom_LSC180(:,3),'Vartype','unequal','Alpha',0.01);
[h365_LSC_AraC25_Mom1_25,p365_LSC_AraC25_Mom1_25] = ttest2(AraC_LSC365(:,2),Mom_LSC365(:,3),'Vartype','unequal','Alpha',0.01);

%HSC
[h90_HSC_AraC125_Bud10,p90_HSC_AraC125_Bud10] = ttest2(AraC_HSC90(:,1),Bud_HSC90(:,2),'Vartype','unequal','Alpha',0.01);
[h180_HSC_AraC125_Bud10,p180_HSC_AraC125_Bud10] = ttest2(AraC_HSC180(:,1),Bud_HSC180(:,2),'Vartype','unequal','Alpha',0.01);
[h365_HSC_AraC125_Bud10,p365_HSC_AraC125_Bud10] = ttest2(AraC_HSC365(:,1),Bud_HSC365(:,2),'Vartype','unequal','Alpha',0.01);

[h90_HSC_AraC125_Mom1_25,p90_HSC_AraC125_Mom1_25] = ttest2(AraC_HSC90(:,1),Mom_HSC90(:,3),'Vartype','unequal','Alpha',0.01);
[h180_HSC_AraC125_Mom1_25,p180_HSC_AraC125_Mom1_25] = ttest2(AraC_HSC180(:,1),Mom_HSC180(:,3),'Vartype','unequal','Alpha',0.01);
[h365_HSC_AraC125_Mom1_25,p365_HSC_AraC125_Mom1_25] = ttest2(AraC_HSC365(:,1),Mom_HSC365(:,3),'Vartype','unequal','Alpha',0.01);

[h90_HSC_AraC25_Bud10,p90_HSC_AraC25_Bud10] = ttest2(AraC_HSC90(:,2),Bud_HSC90(:,2),'Vartype','unequal','Alpha',0.01);
[h180_HSC_AraC25_Bud10,p180_HSC_AraC25_Bud10] = ttest2(AraC_HSC180(:,2),Bud_HSC180(:,2),'Vartype','unequal','Alpha',0.01);
[h365_HSC_AraC25_Bud10,p365_HSC_AraC25_Bud10] = ttest2(AraC_HSC365(:,2),Bud_HSC365(:,2),'Vartype','unequal','Alpha',0.01);

[h90_HSC_AraC25_Mom1_25,p90_HSC_AraC25_Mom1_25] = ttest2(AraC_HSC90(:,2),Mom_HSC90(:,3),'Vartype','unequal','Alpha',0.01);
[h180_HSC_AraC25_Mom1_25,p180_HSC_AraC25_Mom1_25] = ttest2(AraC_HSC180(:,2),Mom_HSC180(:,3),'Vartype','unequal','Alpha',0.01);
[h365_HSC_AraC25_Mom1_25,p365_HSC_AraC25_Mom1_25] = ttest2(AraC_HSC365(:,2),Mom_HSC365(:,3),'Vartype','unequal','Alpha',0.01);

