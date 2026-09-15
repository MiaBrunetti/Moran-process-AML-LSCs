function [fit_parameters,residual,jacobian] = fitting_toxicity(drug_conc,cell_viab,param_guess)

options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','off','MaxFunEval',1000); % setting the optimisation routine specifics

[fit_parameters,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(@residualsfunction, param_guess, [95 0 0], [105 5 10], options); %Invoking optimiser

%------------------------------------------------------------------------
function val = residualsfunction(param)
    
    Emax = param(1);
    IC50 = param(2);
    h    = param(3);
    
    IC50_curve = Emax*((drug_conc.^h)./(IC50.^h+drug_conc.^h));
    %IC50_curve =(Emax-Emin)*IC50./(IC50+drug_conc)+Emin;
    val = IC50_curve - cell_viab;
end
%------------------------------------------------------------------------
end