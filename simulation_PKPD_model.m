function sol = simulation_PKPD_model(p)
%Simulates the pharmacokinetics of the treatment model given the drug, the
%dose, and the frequency of its administration

initialCondition = p.initial_conds;

%Dosing statement
if p.NumberAdmins == 1
    
    p.CurrentAdmin = p.AdministrationTimes(1);
    %Call the solver to calculate solution over initial administration
    sol = ode15s(@(t,A) fitting_PK(A,p), p.TotalTime, initialCondition);

end
if p.NumberAdmins > 1
    
    p.CurrentAdmin = p.AdministrationTimes(1);
    %Call the solver to calculate solution over initial administration
    sol = ode15s(@(t,A) fitting_PK(A,p),[p.AdministrationTimes(1) p.AdministrationTimes(2)], initialCondition);
    
    if numel(p.AdministrationTimes) > 2
        %Extend the solution for each subsequent administration
        for nn = 3:1:numel(p.AdministrationTimes)
            p.CurrentAdmin = p.AdministrationTimes(nn-1);
            %Set ICs to be final solution in last interval
            initialCondition = [p.Dose+sol.y(1,end) sol.y(2:end,end)'];
            sol = odextend(sol,@(t,A) fitting_PK(A,p),p.AdministrationTimes(nn),initialCondition);
        end
        %Calculate solution from end of last admin to final time
        %Set ICs to be final solution in last interval
        p.CurrentAdmin = p.AdministrationTimes(end);
        initialCondition = [p.Dose+sol.y(1,end) sol.y(2:end,end)'];
        sol = odextend(sol,@(t,A) fitting_PK(A,p),p.TotalTime(end),initialCondition);
    else
        %Set ICs to be final solution in last interval
        initialCondition = [p.Dose+sol.y(1,end) sol.y(2:end,end)'];
        p.CurrentAdmin = p.AdministrationTimes(2);
        %Calculate solution from end of last admin to final time
        sol = odextend(sol,@(t,A) fitting_PK(A,p),p.TotalTime(end),initialCondition);
    end

end

end

%------------------------------------------------------------------------
function [eqn] = fitting_PK(A,p)
%Describes a 2 or 3 compartments PK model that has an iv bolus administration.

if p.compartmentNb == 2
    dAc = p.k21_pop*A(2) - (p.k12_pop + p.k_pop)*A(1);
    dA2 = - p.k21_pop*A(2) + p.k12_pop*A(1);

    eqn = [dAc;dA2];
end
    
if p.compartmentNb == 3    
    dAc = p.k21_pop*A(2) + p.k31_pop*A(3) - (p.k12_pop + p.k13_pop + p.k_pop)*A(1);
    dA2 = - p.k21_pop*A(2) + p.k12_pop*A(1);
    dA3 = - p.k31_pop*A(3) + p.k13_pop*A(1);
    
    eqn = [dAc;dA2;dA3];
end

end
%------------------------------------------------------------------------