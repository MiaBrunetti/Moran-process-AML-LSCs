function [eqn] = fitting_ProA_PK(t,A,p)
%The administration is via an infusion (In).
%The PK model has a central compartment (volume V), a peripheral compartment 
%(rate of transfer to and from k12 and k21), and a linear elimination (elimination rate k).

In = interp1(p.PK_simTime,p.In,t);

dAc = In  + p.k21_pop*A(2) - (p.k12_pop + p.k_pop)*A(1);
dA2 = - p.k21_pop*A(2) + p.k12_pop*A(1);

eqn = [dAc;dA2];

end

