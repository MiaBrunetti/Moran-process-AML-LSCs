%% Commands for treatment model

%load work from Commands_PK_fitting.mat, Commands_Cell_Viability_fitting.mat, and Commands_Toxicity_fitting.mat
ProA_PK = load('ProA_PK.mat');
Dig_PK = load('Digoxin_PK.mat');
Oua_PK = load('Ouabain_PK.mat');
Bud_PK = load('Budesonide_PK.mat');
Mom_PK = load('Mometasone_PK.mat');

ProA_cell_viab = load('ProA_cell_viability.mat');
Dig_cell_viab = load('Digoxin_cell_viability.mat');
Oua_cell_viab = load('Ouabain_cell_viability.mat');
Bud_cell_viab = load('Budesonide_cell_viability.mat');
Mom_cell_viab = load('Mometasone_cell_viability.mat');

ProA_tox = load('ProA_toxicity.mat');
Dig_tox = load('Digoxin_toxicity.mat');
Oua_tox = load('Ouabain_toxicity.mat');

%merge structures
ProA = ProA_PK;
f = fieldnames(ProA_cell_viab);
for i = 1:length(f)
    ProA.(f{i}) = ProA_cell_viab.(f{i});
end
g = fieldnames(ProA_tox);
for i = 1:length(g)
    ProA.(g{i}) = ProA_tox.(g{i});
end

Dig = Dig_PK;
f = fieldnames(Dig_cell_viab);
for i = 1:length(f)
    Dig.(f{i}) = Dig_cell_viab.(f{i});
end
g = fieldnames(Dig_tox);
for i = 1:length(g)
    Dig.(g{i}) = Dig_tox.(g{i});
end

Oua = Oua_PK;
f = fieldnames(Oua_cell_viab);
for i = 1:length(f)
    Oua.(f{i}) = Oua_cell_viab.(f{i});
end
g = fieldnames(Oua_tox);
for i = 1:length(g)
    Oua.(g{i}) = Oua_tox.(g{i});
end

Bud = Bud_PK;
f = fieldnames(Bud_cell_viab);
for i = 1:length(f)
    Bud.(f{i}) = Bud_cell_viab.(f{i});
end

Mom = Mom_PK;
f = fieldnames(Mom_cell_viab);
for i = 1:length(f)
    Mom.(f{i}) = Mom_cell_viab.(f{i});
end

%save work
% save('ProA_PKPD.mat', '-struct', 'ProA');
% save('Digoxin_PKPD.mat', '-struct', 'Dig');
% save('Ouabain_PKPD.mat', '-struct', 'Oua');
% save('Budesonide_PKPD.mat', '-struct', 'Bud');
% save('Mometasone_PKPD.mat', '-struct', 'Mom');
