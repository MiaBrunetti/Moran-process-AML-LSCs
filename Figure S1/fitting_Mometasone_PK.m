function [eqn] = fitting_Mometasone_PK(t,A,p)
%The administration is via an infusion (In).
%The PK model has a central compartment (volume V), two peripheral compartments 
%(rate of transfer to and from k12, k21, k13, k31), and a linear elimination (elimination rate k).

%In = interp1(p.PK_simTime,p.In,t);

dAc = p.k21_pop*A(2) +p.k31_pop*A(3) - (p.k12_pop + p.k13_pop + p.k_pop)*A(1);% In  + p.k21_pop*A(2) +p.k31_pop*A(3) - (p.k12_pop + p.k13_pop + p.k_pop)*A(1);
dA2 = - p.k21_pop*A(2) + p.k12_pop*A(1);
dA3 = - p.k31_pop*A(3) + p.k13_pop*A(1);

eqn = [dAc;dA2;dA3];

end

