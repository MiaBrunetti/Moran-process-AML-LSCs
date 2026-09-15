%% Commands for Moran model of stem cells
% Prints all figures for the AraC dynamics and perform t-test.

%load work form Commands_Treatment_Model.mat
% ProA = load('ProA_Moran_Final.mat');
% Dig = load('Dig_Moran_Final.mat');
% Oua = load('Oua_Moran_Final.mat');
% Bud = load('Bud_Moran_Final.mat');
% Mom = load('Mom_Moran_Final.mat');
% AraC = load('AraC_Moran_Final.mat');

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

%% Figure: Cell Viability with AraC

AraC1 = figure;

hold on 
h(1) = scatter(AraC.v_dataConc,AraC.vHSC_data,'o','MarkerEdgeColor','k','MarkerFaceColor','#8EB1DC');
h(2) = errorbar(AraC.v_meanConc,AraC.vHSC_mean,AraC.vHSC_dataLB,AraC.vHSC_dataUB,'^','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',1.0,'Color','b'); %plotting the data
h(3) = plot(AraC.vHSC_simConc,AraC.vHSC_fit,'b','LineWidth',1.5); % plotting curve
h(4) = patch([AraC.vHSC_simConc,fliplr(AraC.vHSC_simConc)],[AraC.vHSCbci_sim(:,1)',fliplr(AraC.vHSCbci_sim(:,2)')],1,'facecolor', 'b', 'edgecolor', 'none', 'facealpha', 0.1); %CI
h(5) = scatter(AraC.v_dataConc,AraC.vLSC_data,'o','MarkerEdgeColor','k','MarkerFaceColor','#E58A8C');
h(6) = errorbar(AraC.v_meanConc,AraC.vLSC_mean,AraC.vLSC_dataLB,AraC.vLSC_dataUB,'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10,'LineWidth',1.0,'Color','r'); %plotting the data
h(7) = plot(AraC.vLSC_simConc,AraC.vLSC_fit,'r','LineWidth',1.5); % plotting curve
h(8) = patch([AraC.vLSC_simConc,fliplr(AraC.vLSC_simConc)],[AraC.vLSCbci_sim(:,1)',fliplr(AraC.vLSCbci_sim(:,2)')],1,'facecolor', 'r', 'edgecolor', 'none', 'facealpha', 0.1); %CI

l(1) = plot(nan,nan,'s','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',10); % plotting curve
l(2) = plot(nan,nan,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',10); % plotting curve
l(3) = plot(nan,nan,'b','LineWidth',1.5); % plotting curve
l(4) = plot(nan,nan,'r','LineWidth',1.5); % plotting curve
hold off
legend(l, {'CD34+ Cord Blood','CD34+ AML 8227','Fit for HSC', 'Fit for LSC'},'Orientation','horizontal','Location', 'southoutside','FontSize',14)
% legend(l, {'CD34- AML 8227','CD34+ Cord Blood','CD34+ AML 8227','Fit for HSC', 'Fit for LSC'},'Orientation','horizontal','Location', 'southoutside','FontSize',14)
% set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
set(gca,'xscale','log','FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
% xlim([0 100])
xlim([1e-02 1000])
ylim([0 105])
xlabel('Concentration (nM)')
ylabel('Cell Viability (%)')
title('AraC','FontSize',20)

%% Figure: Fitness Advantages (AraC)

c_AraC = char('#0B84A5','#F6C85F');
c_AraC = hex2rgb(c_AraC);

AraC2 = figure;

hold on
legend_plot1 = plot(AraC.vHSC_simConc,max(0,AraC.sd),'LineWidth',2,'Color',c_AraC(2,:));
legend_plot2 = plot(AraC.vHSC_simConc,AraC.sp,'LineWidth',2,'Color',c_AraC(1,:));
hold off
legend([legend_plot1 legend_plot2], {'s_d','s_p'},'Orientation','vertical','Location','eastoutside','FontSize',14)
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
set(gca, 'ytick', 0:0.25:1.5);
xlim([0 1000])
ylim([-0.1 1.50])
xlabel('Concentration (nM)')
ylabel('Fitness Value')
set(gca,'xscale','log','FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('AraC','FontSize',20)

%% Figure: Predicted PKPD responses of AraC

AraC3 = figure;
c_AraC = char('#B33DC6','#27AEEF','#87BC45','#EF9B20');
c_AraC = hex2rgb(c_AraC);

hold on
for i = 1:size(AraC.doses,2)
    plot(timeofDiv/60,AraC.TreatmentCc(i,:),'LineWidth',2,'Color',c_AraC(i,:))
end
hold off
legend(sprintf('%g nM',AraC_doses_nM(1,1)),sprintf('%g nM',AraC_doses_nM(1,2)),sprintf('%g nM',AraC_doses_nM(1,3)),sprintf('%g nM',AraC_doses_nM(1,4)),'Location', 'eastoutside','FontSize',16)
xlabel('Time (hours)')
ylabel('Plasma Concentration (nM)')
xlim([0,72])
ylim([1e-05,2e+02])
set(gca,'yscale','log','FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('AraC','FontSize',20)

%% Figure for LSC expansion with AraC

c_AraC = char('#B33DC6','#27AEEF','#87BC45','#EF9B20','#000000');
color_AraC = hex2rgb(c_AraC);
options.legend = [AraC_doses_nM 0];
options.color_area = color_AraC;
options.color_line = color_AraC;
options.alpha = 0.1;
options.line_width = 2;
options.error = 'std';
options.x_axis = timeofDiv/(60*24);

AraC4 = figure;
tiledlayout(2,1,'TileSpacing','compact');

nexttile
plot_areaerrorbar_multiple(AraC.LSC,options)
legend('Orientation','vertical','Location', 'eastoutside','FontSize',18)
ylabel('LSC Number')
xlim([0 sim_time/(60*24)])
ylim([0 35000])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('Cytarabine','FontSize',20)

nexttile
plot_areaerrorbar_multiple(AraC.vLSC,options)
xlabel('Time (days)')
ylabel('LSC Viability (%)')
xlim([0 sim_time/(60*24)])
ylim([0 120])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])

%% Figure for LSC expansion with Budesonide

c_Glu = char('#FF194D','#FFCC0D','#6CADA1','#2A5F65','#000000');
color_Glu = hex2rgb(c_Glu);
options.legend = [Glu_doses_nM 0];
options.color_area = color_Glu;
options.color_line = color_Glu;
options.alpha = 0.1;
options.line_width = 2;
options.error = 'std';
options.x_axis = timeofDiv/(60*24);

fig_Bud = figure;
tiledlayout(2,1,'TileSpacing','compact');

nexttile
plot_areaerrorbar_multiple(Bud.LSC,options)
legend('Orientation','vertical','Location', 'eastoutside','FontSize',18)
ylabel('LSC Number')
xlim([0 sim_time/(60*24)])
ylim([0 35000])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('Budesonide','FontSize',20)

nexttile
plot_areaerrorbar_multiple(Bud.vLSC,options)
xlabel('Time (days)')
ylabel('LSC Viability (%)')
xlim([0 sim_time/(60*24)])
ylim([0 120])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])

%% Figure for LSC expansion with Mometasone

options.legend = [Glu_doses_nM 0];
options.color_area = color_Glu;
options.color_line = color_Glu;
options.alpha = 0.1;
options.line_width = 2;
options.error = 'std';
options.x_axis = timeofDiv/(60*24);

fig_Mom = figure;
tiledlayout(2,1,'TileSpacing','compact');

nexttile
plot_areaerrorbar_multiple(Mom.LSC,options)
legend('Orientation','vertical','Location', 'eastoutside','FontSize',18)
ylabel('LSC Number')
xlim([0 sim_time/(60*24)])
ylim([0 35000])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('Mometasone','FontSize',20)

nexttile
plot_areaerrorbar_multiple(Mom.vLSC,options)
xlabel('Time (days)')
ylabel('LSC Viability (%)')
xlim([0 sim_time/(60*24)])
ylim([0 120])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])

%% Mean value and standard deviation of Moran process at day 90 & 365 (AraC)

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

%% T-test

h_AraC_365 = zeros(1,size(AraC_Moran,1)-1);
p_AraC_365 = zeros(1,size(AraC_Moran,1)-1);
h_AraC_90 = zeros(1,size(AraC_Moran,1)-1);
p_AraC_90 = zeros(1,size(AraC_Moran,1)-1);


for i = 1:size(AraC_Moran,1)-1
    [h_365,p_365] = ttest2(AraC_HSC365(:,5),AraC_HSC365(:,i),'Vartype','unequal','Alpha',0.01);
    h_AraC_365(1,i) = h_365; p_AraC_365(1,i) = p_365;
    [h_90,p_90] = ttest2(AraC_HSC90(:,5),AraC_HSC90(:,i),'Vartype','unequal','Alpha',0.01);
    h_AraC_90(1,i) = h_90; p_AraC_90(1,i) = p_90;
end

%% Boxchart for LSC# at day 90 & 365 (AraC)

box_color = flipud(color_AraC);

AraC5 = figure;
tiledlayout(2,1,'TileSpacing','compact');

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
ylim([1.5e+04 4.5e+04])
ylabel('LSC Number at Day 365')
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
title('Cytarabine','FontSize',20)

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
xlabel('Concentration (nM)')
ylabel('LSC Number at Day 90')
ylim([5000 11000])
set(gca,'FontSize',18,'TickDir','out','TickLength',[0.02 0.025])
ax = gca;
ax.YAxis.Exponent = 3;

%% Code to save figures

% set(AraC1,'Units','Inches');
% pos = get(AraC1,'Position');
% set(AraC1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(AraC1,'Cell Viability AraC','-dpdf','-r0')
% 
% set(AraC2,'Units','Inches');
% pos = get(AraC2,'Position');
% set(AraC2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(AraC2,'Fitness Advantage AraC','-dpdf','-r0')
% 
% set(AraC3,'Units','Inches');
% pos = get(AraC3,'Position');
% set(AraC3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(AraC3,'PK Treatment AraC','-dpdf','-r0')

% set(AraC4,'Units','Inches');
% pos = get(AraC4,'Position');
% set(AraC4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(AraC4,'Moran AraC','-dpdf','-r0')
% 
set(AraC5,'Units','Inches');
pos = get(AraC5,'Position');
set(AraC5,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(AraC5,'T-test AraC','-dpdf','-r0')

%% Comparaison other drugs to cytarabine response

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

%% Save figure

set(bh_HSC125,'Units','Inches');
pos = get(bh_HSC125,'Position');
set(bh_HSC125,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(bh_HSC125,'HSC diff 125nM AraC','-dpdf','-r0')

set(bh_HSC25,'Units','Inches');
pos = get(bh_HSC25,'Position');
set(bh_HSC25,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(bh_HSC25,'HSC diff 25nM AraC','-dpdf','-r0')

set(bh_LSC125,'Units','Inches');
pos = get(bh_LSC125,'Position');
set(bh_LSC125,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(bh_LSC125,'LSC diff 125nM AraC','-dpdf','-r0')

set(bh_LSC25,'Units','Inches');
pos = get(bh_LSC25,'Position');
set(bh_LSC25,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(bh_LSC25,'LSC diff 25nM AraC','-dpdf','-r0')

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


