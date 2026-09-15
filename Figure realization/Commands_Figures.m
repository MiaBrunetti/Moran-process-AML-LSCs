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

%% Figure for LSC expansion and viability

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

%% Figure for LSC expansion with glucocorticoid treatment

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

%% Code to save figures
% 
% set(F2A,'Units','Inches');
% pos = get(F2A,'Position');
% set(F2A,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(F2A,'Cell Viability','-dpdf','-r0')
% 
% set(F2B,'Units','Inches');
% pos = get(F2B,'Position');
% set(F2B,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(F2B,'Fitness Advantage','-dpdf','-r0')
% 
% set(F4A,'Units','Inches');
% pos = get(F4A,'Position');
% set(F4A,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(F4A,'PKPD CarGly','-dpdf','-r0')
% 
% set(F4B,'Units','Inches');
% pos = get(F4B,'Position');
% set(F4B,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(F4B,'PK GluAraC','-dpdf','-r0')
% 
% set(F5A,'Units','Inches');
% pos = get(F5A,'Position');
% set(F5A,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(F5A,'Moran CarGly','-dpdf','-r0')
% 
% set(F5B,'Units','Inches');
% pos = get(F5B,'Position');
% set(F5B,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(F5B,'Moran GluAraC','-dpdf','-r0')

%% Mean value and standard deviation of Moran process at day 90 & 365

ProA_Moran = ProA.LSC;
Dig_Moran = Dig.LSC;
Oua_Moran = Oua.LSC;
Bud_Moran = Bud.LSC;
Mom_Moran = Mom.LSC;

ProA_Via = ProA.vLSC;
Dig_Via = Dig.vLSC;
Oua_Via = Oua.vLSC;
Bud_Via = Bud.vLSC;
Mom_Via = Mom.vLSC;

index365 = 365*24*60*div;
index90 = 90*24*60*div;

% Proscillaridin A
ProA_365 = zeros(100,size(ProA_Moran,1));
ProA_90 = zeros(100,size(ProA_Moran,1));
ProA_v365 = zeros(100,size(ProA_Via,1));
ProA_v90 = zeros(100,size(ProA_Via,1));

for i = 1:size(ProA_Moran,1) %last index is no treatment
    data_LSC = ProA_Moran{i,1};

    % Finding LSC distributions at day 365 and 90
    ProA_365(:,i) = data_LSC(:,index365);
    ProA_90(:,i) = data_LSC(:,index90);
end

for i = 1:size(ProA_Via,1)
    data_vLSC = ProA_Via{i,1};

    % Finding LSC Viability distributions at day 365 and 90
    ProA_v365(:,i) = data_vLSC(:,index365);
    ProA_v90(:,i) = data_vLSC(:,index90);
end

% Digoxin
Dig_365 = zeros(100,size(Dig_Moran,1));
Dig_90 = zeros(100,size(Dig_Moran,1));
Dig_v365 = zeros(100,size(Dig_Via,1));
Dig_v90 = zeros(100,size(Dig_Via,1));

for i = 1:size(Dig_Moran,1) %last index is no treatment
    data_LSC = Dig_Moran{i,1};

    % Finding LSC distributions at day 365 and 90
    Dig_365(:,i) = data_LSC(:,index365);
    Dig_90(:,i) = data_LSC(:,index90);
end

for i = 1:size(Dig_Via,1)
    data_vLSC = Dig_Via{i,1};

    % Finding LSC Viability distributions at day 365 and 90
    Dig_v365(:,i) = data_vLSC(:,index365);
    Dig_v90(:,i) = data_vLSC(:,index90);
end

% Ouabain
Oua_365 = zeros(100,size(Oua_Moran,1));
Oua_90 = zeros(100,size(Oua_Moran,1));
Oua_v365 = zeros(100,size(Oua_Via,1));
Oua_v90 = zeros(100,size(Oua_Via,1));

for i = 1:size(Oua_Moran,1) %last index is no treatment
    data_LSC = Oua_Moran{i,1};

    % Finding LSC distributions at day 365 and 90
    Oua_365(:,i) = data_LSC(:,index365);
    Oua_90(:,i) = data_LSC(:,index90);
end

for i = 1:size(Oua_Via,1)
    data_vLSC = Oua_Via{i,1};

    % Finding LSC Viability distributions at day 365 and 90
    Oua_v365(:,i) = data_vLSC(:,index365);
    Oua_v90(:,i) = data_vLSC(:,index90);
end

% Budesonide
Bud_365 = zeros(100,size(Bud_Moran,1));
Bud_90 = zeros(100,size(Bud_Moran,1));
Bud_v365 = zeros(100,size(Bud_Via,1));
Bud_v90 = zeros(100,size(Bud_Via,1));

for i = 1:size(Bud_Moran,1) %last index is no treatment
    data_LSC = Bud_Moran{i,1};

    % Finding LSC distributions at day 365 and 90
    Bud_365(:,i) = data_LSC(:,index365);
    Bud_90(:,i) = data_LSC(:,index90);
end

for i = 1:size(Bud_Via,1)
    data_vLSC = Bud_Via{i,1};

    % Finding LSC Viability distributions at day 365 and 90
    Bud_v365(:,i) = data_vLSC(:,index365);
    Bud_v90(:,i) = data_vLSC(:,index90);
end

% Mometasone
Mom_365 = zeros(100,size(Mom_Moran,1));
Mom_90 = zeros(100,size(Mom_Moran,1));
Mom_v365 = zeros(100,size(Mom_Via,1));
Mom_v90 = zeros(100,size(Mom_Via,1));

for i = 1:size(Mom_Moran,1) %last index is no treatment
    data_LSC = Mom_Moran{i,1};

    % Finding LSC distributions at day 365 and 90
    Mom_365(:,i) = data_LSC(:,index365);
    Mom_90(:,i) = data_LSC(:,index90);
end

for i = 1:size(Mom_Via,1)
    data_vLSC = Mom_Via{i,1};

    % Finding LSC Viability distributions at day 365 and 90
    Mom_v365(:,i) = data_vLSC(:,index365);
    Mom_v90(:,i) = data_vLSC(:,index90);
end

AraC_Moran = AraC.LSC;
AraC_Via = AraC.vLSC;

index365 = 365*24*60*div;
index90 = 90*24*60*div;

AraC_HSC365 = zeros(100,size(AraC_Moran,1));
AraC_HSC90 = zeros(100,size(AraC_Moran,1));
AraC_v365 = zeros(100,size(AraC_Via,1));
AraC_v90 = zeros(100,size(AraC_Via,1));

for i = 1:size(AraC_Moran,1) %last index is no treatment
    data_LSC = AraC_Moran{i,1};

    % Finding LSC distributions at day 365 and 90
    AraC_HSC365(:,i) = data_LSC(:,index365);
    AraC_HSC90(:,i) = data_LSC(:,index90);
end

for i = 1:size(AraC_Via,1)
    data_vLSC = AraC_Via{i,1};

    % Finding LSC Viability distributions at day 365 and 90
    AraC_v365(:,i) = data_vLSC(:,index365);
    AraC_v90(:,i) = data_vLSC(:,index90);
end

%% Boxchart for LSC# and viability at day 90 & 365 with cardiac glycoside treatment

box_colorCarGly = flipud(color_CarGly);

b_CarGly = figure;
tiledlayout(2,3,'TileSpacing','compact');

nexttile %Proscillaridin A (day 90)
fProA_90 = fliplr(ProA_90);
[N,M] = size(fProA_90);
hold on
for i = 1:M
    boxchart(fProA_90(:,i),'MarkerStyle','none','BoxFaceColor',box_colorCarGly(i,:),'XData',i*ones(N,1))
    plot(i,mean(fProA_90(:,i)),'*','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor',box_colorCarGly(i,:))
end
ax = gca();
ax.XAxis.Categories = categorical(1:M);
xticklabels({'0','10','20','30','50'})
ylabel('LSC Number at Day 90')
ylim([4500 10000])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
ax = gca;
ax.YAxis.Exponent = 3;
title('Proscillaridin A','FontSize',20)

nexttile %Digoxin (day 90)
fDig_90 = fliplr(Dig_90);
[N,M] = size(fDig_90);
hold on
for i = 1:M
    boxchart(fDig_90(:,i),'MarkerStyle','none','BoxFaceColor',box_colorCarGly(i,:),'XData',i*ones(N,1))
    plot(i,mean(fDig_90(:,i)),'*','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor',box_colorCarGly(i,:))
end
ax = gca();
ax.XAxis.Categories = categorical(1:M);
xticklabels({'0','10','20','30','50'})
ylim([3750 10000])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
ax = gca;
ax.YAxis.Exponent = 3;
title('Digoxin','FontSize',20)

nexttile %Ouabain (day 90)
fOua_90 = fliplr(Oua_90);
[N,M] = size(fOua_90);
hold on
for i = 1:M
    boxchart(fOua_90(:,i),'MarkerStyle','none','BoxFaceColor',box_colorCarGly(i,:),'XData',i*ones(N,1))
    plot(i,mean(fOua_90(:,i)),'*','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor',box_colorCarGly(i,:))
end
ax = gca();
ax.XAxis.Categories = categorical(1:M);
xticklabels({'0','10','20','30','50'})
ylim([6250 10000])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
ax = gca;
ax.YAxis.Exponent = 3;
title('Ouabain','FontSize',20)

nexttile %Proscillaridin A (day 365)
fProA_365 = fliplr(ProA_365);
[N,M] = size(fProA_365);
hold on
for i = 1:M
    boxchart(fProA_365(:,i),'MarkerStyle','none','BoxFaceColor',box_colorCarGly(i,:),'XData',i*ones(N,1))
    plot(i,mean(fProA_365(:,i)),'*','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor',box_colorCarGly(i,:))
end
ax = gca();
ax.XAxis.Categories = categorical(1:M);
xticklabels({'0','10','20','30','50'})
ylabel('LSC Number at Day 365')
ylim([9750 4e+04])
xlabel('Concentration (nM)')
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])

nexttile %Digoxin (day 365)
fDig_365 = fliplr(Dig_365);
[N,M] = size(fDig_365);
hold on
for i = 1:M
    boxchart(fDig_365(:,i),'MarkerStyle','none','BoxFaceColor',box_colorCarGly(i,:),'XData',i*ones(N,1))
        plot(i,mean(fDig_365(:,i)),'*','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor',box_colorCarGly(i,:))
end
ax = gca();
ax.XAxis.Categories = categorical(1:M);
xticklabels({'0','10','20','30','50'})
ylim([0.4e+04 4e+04])
xlabel('Concentration (nM)')
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])

nexttile %Ouabain (day 365)
fOua_365 = fliplr(Oua_365);
[N,M] = size(Oua_365);
hold on
for i = 1:M
    boxchart(fOua_365(:,i),'MarkerStyle','none','BoxFaceColor',box_colorCarGly(i,:),'XData',i*ones(N,1))
        plot(i,mean(fOua_365(:,i)),'*','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor',box_colorCarGly(i,:))
end
ax = gca();
ax.XAxis.Categories = categorical(1:M);
xticklabels({'0','10','20','30','50'})
ylim([2.4e+04 4e+04])
xlabel('Concentration (nM)')
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])

%% Boxchart for LSC# and viability at day 90 & 365 with glucocorticoid treatment

c_AraC = char('#B33DC6','#27AEEF','#87BC45','#EF9B20','#000000');
color_AraC = hex2rgb(c_AraC);

box_colorGlu = flipud(color_Glu);
box_color = flipud(color_AraC);

b_Glu = figure;
tiledlayout(2,3,'TileSpacing','compact');

nexttile %Budesonide (day 90)
fBud_90 = fliplr(Bud_90);
[N,M] = size(fBud_90);
hold on
for i = 1:M
    boxchart(fBud_90(:,i),'MarkerStyle','none','BoxFaceColor',box_colorGlu(i,:),'XData',i*ones(N,1))
    plot(i,mean(fBud_90(:,i)),'*','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor',box_colorGlu(i,:))
end
ax = gca();
ax.XAxis.Categories = categorical(1:M);
xticklabels({'0','0.25','1.5','10','25'})
ylabel('LSC Number at Day 90')
ylim([3750 11000])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
ax = gca;
ax.YAxis.Exponent = 3;
title('Budesonide','FontSize',20)

nexttile %Mometasone (day 90)
fMom_90 = fliplr(Mom_90);
[N,M] = size(fMom_90);
hold on
for i = 1:M
    boxchart(fMom_90(:,i),'MarkerStyle','none','BoxFaceColor',box_colorGlu(i,:),'XData',i*ones(N,1))
    plot(i,mean(fMom_90(:,i)),'*','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor',box_colorGlu(i,:))
end
ax = gca();
ax.XAxis.Categories = categorical(1:M);
xticklabels({'0','0.25','1.5','10','25'})
ylim([-50 11000])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
ax = gca;
ax.YAxis.Exponent = 3;
title('Mometasone','FontSize',20)

nexttile %AraC (day 90)
fAraC_90 = fliplr(AraC_HSC90);
[N,M] = size(fAraC_90);
hold on
for i = 1:M
    boxchart(fAraC_90(:,i),'MarkerStyle','none','BoxFaceColor',box_color(i,:),'XData',i*ones(N,1))
    plot(i,mean(fAraC_90(:,i)),'*','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor',box_color(i,:))
end
ax = gca();
ax.XAxis.Categories = categorical(1:M);
xticklabels({'0','1','5','25','125'})
ylabel('LSC Number at Day 90')
ylim([5000 11000])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
ax = gca;
ax.YAxis.Exponent = 3;
title('Cytarabine','FontSize',20)

nexttile %Budesonide (day 365)
fBud_365 = fliplr(Bud_365);
[N,M] = size(fBud_365);
hold on
for i = 1:M
    boxchart(fBud_365(:,i),'MarkerStyle','none','BoxFaceColor',box_colorGlu(i,:),'XData',i*ones(N,1))
    plot(i,mean(fBud_365(:,i)),'*','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor',box_colorGlu(i,:))
end
ax = gca();
ax.XAxis.Categories = categorical(1:M);
xticklabels({'0','0.25','1.5','10','25'})
xlabel('Concentration (nM)')
ylabel('LSC Number at Day 365')
ylim([3500 4.5e+04])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])

nexttile %Mometasone (day 365)
fMom_365 = fliplr(Mom_365);
[N,M] = size(fMom_365);
hold on
for i = 1:M
    boxchart(fMom_365(:,i),'MarkerStyle','none','BoxFaceColor',box_colorGlu(i,:),'XData',i*ones(N,1))
    plot(i,mean(fMom_365(:,i)),'*','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor',box_colorGlu(i,:))
end
ax = gca();
ax.XAxis.Categories = categorical(1:M);
xticklabels({'0','0.25','1.5','10','25'})
xlabel('Concentration (nM)')
ylim([-500 4.5e+04])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])

nexttile %AraC (day 365)
fAraC_365 = fliplr(AraC_HSC365);
[N,M] = size(fAraC_365);
hold on
for i = 1:M
    boxchart(fAraC_365(:,i),'BoxFaceColor',box_color(i,:),'XData',i*ones(N,1))
    plot(i,mean(fAraC_365(:,i)),'*','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor',box_color(i,:))
end
ax = gca();
ax.XAxis.Categories = categorical(1:M);
xticklabels({'0','1','5','25','125'})
xlabel('Concentration (nM)')
ylim([1.5e+04 4.5e+04])
ylabel('LSC Number at Day 365')
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])

%% Save

set(b_CarGly,'Units','Inches');
pos = get(b_CarGly,'Position');
set(b_CarGly,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(b_CarGly,'T-test CarGly','-dpdf','-r0')

set(b_Glu,'Units','Inches');
pos = get(b_Glu,'Position');
set(b_Glu,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(b_Glu,'T-test Glu AraC','-dpdf','-r0')

