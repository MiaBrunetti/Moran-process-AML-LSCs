%% Commands for Moran model of stem cells

%load work form Commands_Treatment_Model.mat
ProA = load('ProA_TreatmentModel.mat');
Dig = load('Digoxin_TreatmentModel.mat');
Oua = load('Ouabain_TreatmentModel.mat');
Bud = load('Budesonide_TreatmentModel.mat');
Mom = load('Mometasone_TreatmentModel.mat');

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
CarGly_doses_nM = [50 30 20 10]; %in nM
Glu_doses_nM = [25 10 1.5 0.25]; %in nM

%% Proscillaridin A

ProA.LSC = cell(size(CarGly_doses_nM,2)+1,1);
ProA.vLSC = cell(size(CarGly_doses_nM,2),1);

%Simulation without drug
ProA.noDrug_vHSC = ProA.vHSC_Emax + ProA.vHSC_Emin;
ProA.noDrug_vLSC = ProA.vLSC_Emax + ProA.vLSC_Emin;
ProA.noDrug_sd = 0;
ProA.noDrug_sp = 1.25;

for sim = 1:sim_num  
    id = frac*Nd;
    x_LSC = id;
    x_HSC = Nd-x_LSC;
    LSC_diff_pool(sim,1) = 0;
    
    for eventcount = 1:maxNoEvents
       LSC_frequency(sim,eventcount) = x_LSC/x_HSC;
       ProA.noDrug_LSC(sim,eventcount) = x_LSC;

       ProA.d_LSC = x_LSC/(x_LSC+(1-ProA.noDrug_sd)*x_HSC);
       ProA.d_HSC = (1-ProA.noDrug_sd)*x_HSC/(x_LSC+(1-ProA.noDrug_sd)*x_HSC);
       ProA.p_LSC = ProA.noDrug_sp*x_LSC/(ProA.noDrug_sp*x_LSC+x_HSC);
       ProA.p_HSC = x_HSC/(ProA.noDrug_sp*x_LSC+x_HSC);

       event1 = ProA.d_LSC*ProA.p_LSC; 
       event2 = ProA.d_LSC*ProA.p_HSC;
       event3 = ProA.d_HSC*ProA.p_LSC;
       event4 = ProA.d_HSC*ProA.p_HSC;

       eventvec = cumsum([event1,event2,event3,event4]);
       pval = rand;

       if pval<eventvec(1)
          id = id;
          LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount) + 1;
       elseif pval<=eventvec(2)
          id = id-1;
          LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount) + 1;   
       elseif pval<=eventvec(3)
          id = id+1;
          LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount);
       else
          id = id;
          LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount);
       end
       x_LSC = id;
       x_HSC = Nd - x_LSC;
    end
end

ProA.LSC{size(CarGly_doses_nM,2)+1,1} = ProA.noDrug_LSC;
disp('ProA No Drug Loop')

%Simulation with treatment: Every 24hrs (once a day)
ProA.Treatment_sd = 1-ProA.Treatment_vLSC./ProA.Treatment_vHSC - (1-ProA.noDrug_vLSC./ProA.noDrug_vHSC);
ProA.Treatment_sp = 1.25*(ProA.Treatment_vLSC/ProA.noDrug_vLSC);

for i = 1:size(CarGly_doses_nM,2)
    for sim = 1:sim_num
        id = frac*Nd;
        x_LSC = id;
        x_HSC = Nd-x_LSC;
        LSC_diff_pool(sim,1) = 0;

        for eventcount = 1:maxNoEvents
           LSC_frequency(sim,eventcount) = x_LSC/x_HSC;
           ProA.Treatment_LSC(sim,eventcount) = x_LSC;

           %Picking the sd and sp value for the specific division event
           ProA.event_sd = ProA.Treatment_sd(i,eventcount);
           ProA.event_sp = ProA.Treatment_sp(i,eventcount);

           ProA.d_LSC = x_LSC/(x_LSC+(1-ProA.event_sd)*x_HSC);
           ProA.d_HSC = (1-ProA.event_sd)*x_HSC/(x_LSC+(1-ProA.event_sd)*x_HSC);
           ProA.p_LSC = ProA.event_sp*x_LSC/(ProA.event_sp*x_LSC+x_HSC);
           ProA.p_HSC = x_HSC/(ProA.event_sp*x_LSC+x_HSC);

           event1 = ProA.d_LSC*ProA.p_LSC; 
           event2 = ProA.d_LSC*ProA.p_HSC;
           event3 = ProA.d_HSC*ProA.p_LSC;
           event4 = ProA.d_HSC*ProA.p_HSC;

           eventvec = cumsum([event1,event2,event3,event4]);
           pval = rand;

           if pval<eventvec(1)
               id = id;
               LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount) + 1;
           elseif pval<=eventvec(2)
               id = id-1;
               LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount) + 1;   
           elseif pval<=eventvec(3)
               id = id+1;
               LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount);
           else
               id = id;
               LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount);
           end
           x_LSC = id;
           x_HSC = Nd - x_LSC;
        end
    end
    ProA.LSC{i,1} = ProA.Treatment_LSC;
    ProA.vLSC{i,1} = (ProA.Treatment_LSC./ProA.noDrug_LSC).*100;
    disp(['ProA Treatment Loop #' num2str(i) ''])
end

%% Digoxin

Dig.LSC = cell(size(CarGly_doses_nM,2)+1,1); %number of LSCs
Dig.vLSC = cell(size(CarGly_doses_nM,2),1);  %LSC viability: #LSCs with drug/#LSCs without drug

%Simulation without drug
Dig.noDrug_vHSC = Dig.vHSC_Emax + Dig.vHSC_Emin;
Dig.noDrug_vLSC = Dig.vLSC_Emax + Dig.vLSC_Emin;
Dig.noDrug_sd = 0;
Dig.noDrug_sp = 1.25;

for sim = 1:sim_num  
    id = frac*Nd;
    x_LSC = id;
    x_HSC = Nd-x_LSC;
    LSC_diff_pool(sim,1) = 0;
    
    for eventcount = 1:maxNoEvents
       LSC_frequency(sim,eventcount) = x_LSC/x_HSC;
       Dig.noDrug_LSC(sim,eventcount) = x_LSC;

       Dig.d_LSC = x_LSC/(x_LSC+(1-Dig.noDrug_sd)*x_HSC);
       Dig.d_HSC = (1-Dig.noDrug_sd)*x_HSC/(x_LSC+(1-Dig.noDrug_sd)*x_HSC);
       Dig.p_LSC = Dig.noDrug_sp*x_LSC/(Dig.noDrug_sp*x_LSC+x_HSC);
       Dig.p_HSC = x_HSC/(Dig.noDrug_sp*x_LSC+x_HSC);

       event1 = Dig.d_LSC*Dig.p_LSC; 
       event2 = Dig.d_LSC*Dig.p_HSC;
       event3 = Dig.d_HSC*Dig.p_LSC;
       event4 = Dig.d_HSC*Dig.p_HSC;

       eventvec = cumsum([event1,event2,event3,event4]);
       pval = rand;

       if pval<eventvec(1)
          id = id;
          LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount) + 1;
       elseif pval<=eventvec(2)
          id = id-1;
          LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount) + 1;   
       elseif pval<=eventvec(3)
          id = id+1;
          LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount);
       else
          id = id;
          LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount);
       end
       x_LSC = id;
       x_HSC = Nd-x_LSC;
    end
end

Dig.LSC{size(CarGly_doses_nM,2)+1,1} = Dig.noDrug_LSC;
disp('Digoxin No Drug Loop')

%Simulation with treatment: Every 24hrs (once a day)
Dig.Treatment_sd = 1-Dig.Treatment_vLSC./Dig.Treatment_vHSC  - (1-Dig.noDrug_vLSC./Dig.noDrug_vHSC);
Dig.Treatment_sp = 1.25*(Dig.Treatment_vLSC/Dig.noDrug_vLSC);

for i = 1:size(CarGly_doses_nM,2)
    for sim = 1:sim_num
        id = frac*Nd;
        x_LSC = id;
        x_HSC = Nd-x_LSC;
        LSC_diff_pool(sim,1) = 0;

        for eventcount = 1:maxNoEvents
           LSC_frequency(sim,eventcount) = x_LSC/x_HSC;
           Dig.Treatment_LSC(sim,eventcount) = x_LSC;

           %Picking the sd and sp value for the specific division event
           Dig.event_sd = Dig.Treatment_sd(i,eventcount);
           Dig.event_sp = Dig.Treatment_sp(i,eventcount);

           Dig.d_LSC = x_LSC/(x_LSC+(1-Dig.event_sd)*x_HSC);
           Dig.d_HSC = (1-Dig.event_sd)*x_HSC/(x_LSC+(1-Dig.event_sd)*x_HSC);
           Dig.p_LSC = Dig.event_sp*x_LSC/(Dig.event_sp*x_LSC+x_HSC);
           Dig.p_HSC = x_HSC/(Dig.event_sp*x_LSC+x_HSC);

           event1 = Dig.d_LSC*Dig.p_LSC; 
           event2 = Dig.d_LSC*Dig.p_HSC;
           event3 = Dig.d_HSC*Dig.p_LSC;
           event4 = Dig.d_HSC*Dig.p_HSC;

           eventvec = cumsum([event1,event2,event3,event4]);
           pval = rand;

           if pval<eventvec(1)
               id = id;
               LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount) + 1;
           elseif pval<=eventvec(2)
               id = id-1;
               LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount) + 1;   
           elseif pval<=eventvec(3)
               id = id+1;
               LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount);
           else
               id = id;
               LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount);
           end
           x_LSC = id;
           x_HSC = Nd-x_LSC;
        end
    end
    Dig.LSC{i,1} = Dig.Treatment_LSC;
    Dig.vLSC{i,1} = (Dig.Treatment_LSC./Dig.noDrug_LSC).*100;
    disp(['Digoxin Treatment Loop #' num2str(i) ''])
end

%% Ouabain

Oua.LSC = cell(size(CarGly_doses_nM,2)+1,1); %number of LSCs
Oua.vLSC = cell(size(CarGly_doses_nM,2),1);  %LSC viability: #LSCs with drug/#LSCs without drug

%Simulation without drug
Oua.noDrug_vHSC = Oua.vHSC_Emax + Oua.vHSC_Emin;
Oua.noDrug_vLSC = Oua.vLSC_Emax + Oua.vLSC_Emin;
Oua.noDrug_sd = 0;
Oua.noDrug_sp = 1.25;

for sim = 1:sim_num  
    id = frac*Nd;
    x_LSC = id;
    x_HSC = Nd-x_LSC;
    LSC_diff_pool(sim,1) = 0;
    
    for eventcount = 1:maxNoEvents
       LSC_frequency(sim,eventcount) = x_LSC/x_HSC;
       Oua.noDrug_LSC(sim,eventcount) = x_LSC;

       Oua.d_LSC = x_LSC/(x_LSC+(1-Oua.noDrug_sd)*x_HSC);
       Oua.d_HSC = (1-Oua.noDrug_sd)*x_HSC/(x_LSC+(1-Oua.noDrug_sd)*x_HSC);
       Oua.p_LSC = Oua.noDrug_sp*x_LSC/(Oua.noDrug_sp*x_LSC+x_HSC);
       Oua.p_HSC = x_HSC/(Oua.noDrug_sp*x_LSC+x_HSC);

       event1 = Oua.d_LSC*Oua.p_LSC; 
       event2 = Oua.d_LSC*Oua.p_HSC;
       event3 = Oua.d_HSC*Oua.p_LSC;
       event4 = Oua.d_HSC*Oua.p_HSC;

       eventvec = cumsum([event1,event2,event3,event4]);
       pval = rand;

       if pval<eventvec(1)
          id = id;
          LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount) + 1;
       elseif pval<=eventvec(2)
          id = id-1;
          LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount) + 1;   
       elseif pval<=eventvec(3) 
          id = id+1;
          LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount);
       else
          id = id;
          LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount);
       end
       x_LSC = id;
       x_HSC = Nd-x_LSC;
    end
end

Oua.LSC{size(CarGly_doses_nM,2)+1,1} = Oua.noDrug_LSC;
disp('Ouabain No Drug Loop')

%Simulation with treatment: Every 24hrs (once a day)
Oua.Treatment_sd = 1-Oua.Treatment_vLSC./Oua.Treatment_vHSC - (1-Oua.noDrug_vLSC./Oua.noDrug_vHSC);
Oua.Treatment_sp = 1.25*(Oua.Treatment_vLSC/Oua.noDrug_vLSC);

for i = 1:size(CarGly_doses_nM,2)
    for sim = 1:sim_num
        id = frac*Nd;
        x_LSC = id;
        x_HSC = Nd-x_LSC;
        LSC_diff_pool(sim,1) = 0;

        for eventcount = 1:maxNoEvents
           LSC_frequency(sim,eventcount) = x_LSC/x_HSC;
           Oua.Treatment_LSC(sim,eventcount) = x_LSC;

           %Picking the sd and sp value for the specific division event
           Oua.event_sd = Oua.Treatment_sd(i,eventcount);
           Oua.event_sp = Oua.Treatment_sp(i,eventcount);

           Oua.d_LSC = x_LSC/(x_LSC+(1-Oua.event_sd)*x_HSC);
           Oua.d_HSC = (1-Oua.event_sd)*x_HSC/(x_LSC+(1-Oua.event_sd)*x_HSC);
           Oua.p_LSC = Oua.event_sp*x_LSC/(Oua.event_sp*x_LSC+x_HSC);
           Oua.p_HSC = x_HSC/(Oua.event_sp*x_LSC+x_HSC);

           event1 = Oua.d_LSC*Oua.p_LSC; 
           event2 = Oua.d_LSC*Oua.p_HSC;
           event3 = Oua.d_HSC*Oua.p_LSC;
           event4 = Oua.d_HSC*Oua.p_HSC;

           eventvec = cumsum([event1,event2,event3,event4]);
           pval = rand;

           if pval<eventvec(1)
               id = id;
               LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount) + 1;
           elseif pval<=eventvec(2)
               id = id-1;
               LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount) + 1;   
           elseif pval<=eventvec(3)
               id = id+1;
               LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount);
           else
               id = id;
               LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount);
           end
           x_LSC = id;
           x_HSC = Nd-x_LSC;
        end
    end
    Oua.LSC{i,1} = Oua.Treatment_LSC;
    Oua.vLSC{i,1} = (Oua.Treatment_LSC./Oua.noDrug_LSC).*100;
    disp(['Ouabain Treatment Loop #' num2str(i) ''])
end

%% Budesonide

Bud.LSC = cell(size(Glu_doses_nM,2)+1,1); %number of LSCs
Bud.vLSC = cell(size(Glu_doses_nM,2),1);  %LSC viability: #LSCs with drug/#LSCs without drug

%Simulation without drug
Bud.noDrug_vHSC = Bud.vHSC_Emax + Bud.vHSC_Emin;
Bud.noDrug_vLSC = Bud.vLSC_Emax + Bud.vLSC_Emin;
Bud.noDrug_sd = 0;
Bud.noDrug_sp = 1.25;

for sim = 1:sim_num  
    id = frac*Nd;
    x_LSC = id;
    x_HSC = Nd-x_LSC;
    LSC_diff_pool(sim,1) = 0;
    
    for eventcount = 1:maxNoEvents
       LSC_frequency(sim,eventcount) = x_LSC/x_HSC;
       Bud.noDrug_LSC(sim,eventcount) = x_LSC;

       Bud.d_LSC = x_LSC/(x_LSC+(1-Bud.noDrug_sd)*x_HSC);
       Bud.d_HSC = (1-Bud.noDrug_sd)*x_HSC/(x_LSC+(1-Bud.noDrug_sd)*x_HSC);
       Bud.p_LSC = Bud.noDrug_sp*x_LSC/(Bud.noDrug_sp*x_LSC+x_HSC);
       Bud.p_HSC = x_HSC/(Bud.noDrug_sp*x_LSC+x_HSC);

       event1 = Bud.d_LSC*Bud.p_LSC; 
       event2 = Bud.d_LSC*Bud.p_HSC;
       event3 = Bud.d_HSC*Bud.p_LSC;
       event4 = Bud.d_HSC*Bud.p_HSC;

       eventvec = cumsum([event1,event2,event3,event4]);
       pval = rand;

       if pval<eventvec(1)
          id = id;
          LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount) + 1;
       elseif pval<=eventvec(2)
          id = id-1;
          LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount) + 1;   
       elseif pval<=eventvec(3)
          id = id+1;
          LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount);
       else
          id = id;
          LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount);
       end
       x_LSC = id;
       x_HSC = Nd-x_LSC;
    end
end

Bud.LSC{size(Glu_doses_nM,2)+1,1} = Bud.noDrug_LSC;
disp('Budesonide No Drug Loop')

%Simulation with treatment: Every 24hrs (once a day)
Bud.Treatment_sd = 1-Bud.Treatment_vLSC./Bud.Treatment_vHSC  - (1-Bud.noDrug_vLSC./Bud.noDrug_vHSC);
Bud.Treatment_sp = 1.25*(Bud.Treatment_vLSC/Bud.noDrug_vLSC);

for i = 1:size(Glu_doses_nM,2)
    for sim = 1:sim_num
        id = frac*Nd;
        x_LSC = id;
        x_HSC = Nd-x_LSC;
        LSC_diff_pool(sim,1) = 0;

        for eventcount = 1:maxNoEvents
           LSC_frequency(sim,eventcount) = x_LSC/x_HSC;
           Bud.Treatment_LSC(sim,eventcount) = x_LSC;

           %Picking the sd and sp value for the specific division event
           Bud.event_sd = Bud.Treatment_sd(i,eventcount);
           Bud.event_sp = Bud.Treatment_sp(i,eventcount);

           Bud.d_LSC = x_LSC/(x_LSC+(1-Bud.event_sd)*x_HSC);
           Bud.d_HSC = (1-Bud.event_sd)*x_HSC/(x_LSC+(1-Bud.event_sd)*x_HSC);
           Bud.p_LSC = Bud.event_sp*x_LSC/(Bud.event_sp*x_LSC+x_HSC);
           Bud.p_HSC = x_HSC/(Bud.event_sp*x_LSC+x_HSC);

           event1 = Bud.d_LSC*Bud.p_LSC; 
           event2 = Bud.d_LSC*Bud.p_HSC;
           event3 = Bud.d_HSC*Bud.p_LSC;
           event4 = Bud.d_HSC*Bud.p_HSC;

           eventvec = cumsum([event1,event2,event3,event4]);
           pval = rand;

           if pval<eventvec(1)
               id = id;
               LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount)+1;
           elseif pval<=eventvec(2)
               id = id-1;
               LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount)+1;   
           elseif pval<=eventvec(3)
               id = id+1;
               LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount);
           else
               id = id;
               LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount);
           end
           x_LSC = id;
           x_HSC = Nd-x_LSC;
        end
    end
    Bud.LSC{i,1} = Bud.Treatment_LSC;
    Bud.vLSC{i,1} = (Bud.Treatment_LSC./Bud.noDrug_LSC).*100;
    disp(['Budesonide Treatment Loop #' num2str(i) ''])
end

%% Mometasone

Mom.LSC = cell(size(Glu_doses_nM,2)+1,1); %number of LSCs
Mom.vLSC = cell(size(Glu_doses_nM,2),1);  %LSC viability: #LSCs with drug/#LSCs without drug

%Simulation without drug
Mom.noDrug_vHSC = Mom.vHSC_Emax + Mom.vHSC_Emin;
Mom.noDrug_vLSC = Mom.vLSC_Emax + Mom.vLSC_Emin;
Mom.noDrug_sd = 0;
Mom.noDrug_sp = 1.25;

for sim = 1:sim_num  
    id = frac*Nd;
    x_LSC = id;
    x_HSC = Nd-x_LSC;
    LSC_diff_pool(sim,1) = 0;
    
    for eventcount = 1:maxNoEvents
       LSC_frequency(sim,eventcount) = x_LSC/x_HSC;
       Mom.noDrug_LSC(sim,eventcount) = x_LSC;

       Mom.d_LSC = x_LSC/(x_LSC+(1-Mom.noDrug_sd)*x_HSC);
       Mom.d_HSC = (1-Mom.noDrug_sd)*x_HSC/(x_LSC+(1-Mom.noDrug_sd)*x_HSC);
       Mom.p_LSC = Mom.noDrug_sp*x_LSC/(Mom.noDrug_sp*x_LSC+x_HSC);
       Mom.p_HSC = x_HSC/(Mom.noDrug_sp*x_LSC+x_HSC);

       event1 = Mom.d_LSC*Mom.p_LSC; 
       event2 = Mom.d_LSC*Mom.p_HSC;
       event3 = Mom.d_HSC*Mom.p_LSC;
       event4 = Mom.d_HSC*Mom.p_HSC;

       eventvec = cumsum([event1,event2,event3,event4]);
       pval = rand;

       if pval<eventvec(1)
          id = id;
          LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount) + 1;
       elseif pval<=eventvec(2)
          id = id-1;
          LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount) + 1;   
       elseif pval<=eventvec(3)
          id = id+1;
          LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount);
       else
          id = id;
          LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount);
       end
       x_LSC = id;
       x_HSC = Nd-x_LSC;
    end
end

Mom.LSC{size(Glu_doses_nM,2)+1,1} = Mom.noDrug_LSC;
disp('Mometasone No Drug Loop')

%Simulation with treatment: Every 24hrs (once a day)
Mom.Treatment_sd = 1-Mom.Treatment_vLSC./Mom.Treatment_vHSC  - (1-Mom.noDrug_vLSC./Mom.noDrug_vHSC);
Mom.Treatment_sp = 1.25*(Mom.Treatment_vLSC/Mom.noDrug_vLSC);

for i = 1:size(Glu_doses_nM,2)
    for sim = 1:sim_num
        id = frac*Nd;
        x_LSC = id;
        x_HSC = Nd-x_LSC;
        LSC_diff_pool(sim,1) = 0;

        for eventcount = 1:maxNoEvents
           LSC_frequency(sim,eventcount) = x_LSC/x_HSC;
           Mom.Treatment_LSC(sim,eventcount) = x_LSC;

           %Picking the sd and sp value for the specific division event
           Mom.event_sd = Mom.Treatment_sd(i,eventcount);
           Mom.event_sp = Mom.Treatment_sp(i,eventcount);

           Mom.d_LSC = x_LSC/(x_LSC+(1-Mom.event_sd)*x_HSC);
           Mom.d_HSC = (1-Mom.event_sd)*x_HSC/(x_LSC+(1-Mom.event_sd)*x_HSC);
           Mom.p_LSC = Mom.event_sp*x_LSC/(Mom.event_sp*x_LSC+x_HSC);
           Mom.p_HSC = x_HSC/(Mom.event_sp*x_LSC+x_HSC);

           event1 = Mom.d_LSC*Mom.p_LSC; 
           event2 = Mom.d_LSC*Mom.p_HSC;
           event3 = Mom.d_HSC*Mom.p_LSC;
           event4 = Mom.d_HSC*Mom.p_HSC;

           eventvec = cumsum([event1,event2,event3,event4]);
           pval = rand;

           if pval<eventvec(1)
               id = id;
               LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount) + 1;
           elseif pval<=eventvec(2)
               id = id-1;
               LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount) + 1;   
           elseif pval<=eventvec(3)
               id = id+1;
               LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount);
           else
               id = id;
               LSC_diff_pool(sim,eventcount+1) = LSC_diff_pool(sim,eventcount);
           end
           x_LSC = id;
           x_HSC = Nd-x_LSC;
        end
    end
    Mom.LSC{i,1} = Mom.Treatment_LSC;
    Mom.vLSC{i,1} = (Mom.Treatment_LSC./Mom.noDrug_LSC).*100;
    disp(['Mometasone Treatment Loop #' num2str(i) ''])
end

%% Save work

% save('ProA_Moran.mat', '-struct', 'ProA');
% save('Dig_Moran.mat', '-struct', 'Dig');
% save('Oua_Moran.mat', '-struct', 'Oua');
% save('Bud_Moran.mat', '-struct', 'Bud');
% save('Mom_Moran.mat', '-struct', 'Mom');

%% Figure for LSC expansion with cardiac glycoside treatment

c = char('#702A8C','#BF2669','#FF7326','#FFCC0D','#000000');
color = hex2rgb(c);
options.color_area = color;
options.color_line = color;
options.legend = [CarGly_doses_nM 0];
options.alpha = 0.1;
options.line_width = 2;
options.error = 'std';
options.x_axis = timeofDiv/(60*24);

figure
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

c = char('#FF194D','#FFCC0D','#6CADA1','#2A5F65','#000000');
color = hex2rgb(c);
options.legend = [Glu_doses_nM 0];
options.color_area = color;
options.color_line = color;
options.alpha = 0.1;
options.line_width = 2;
options.error = 'std';
options.x_axis = timeofDiv/(60*24);

figure
tiledlayout(2,2,'TileSpacing','compact');

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
