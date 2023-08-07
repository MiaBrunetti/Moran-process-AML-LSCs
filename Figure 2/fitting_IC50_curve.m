function [fit_parameters,residual,jacobian] = fitting_IC50_curve(drug_conc,cell_viab,param_guess)

options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','off','MaxFunEval',1000); % setting the optimisation routine specifics

[fit_parameters,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(@residualsfunction, param_guess, [0 95 0 0], [1 105 100 25], options);   % Invoking optimiser
%ci = nlparci(fit_parameters,residual,'jacobian',jacobian);

%------------------------------------------------------------------------
function val = residualsfunction(param)
    
    Emin = param(1);
    Emax = param(2);
    IC50 = param(3);
    h    = param(4);
    
    IC50_curve = Emax - (Emax*drug_conc.^h)./(IC50.^h+drug_conc.^h) + Emin;
    %IC50_curve =(Emax-Emin)*IC50./(IC50+drug_conc)+Emin;
    val = IC50_curve - cell_viab;
end
%------------------------------------------------------------------------
end