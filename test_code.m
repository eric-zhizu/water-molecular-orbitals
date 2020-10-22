%% Comparison to Illinois state data
% My data
my_energies = [-19.648070410709670;-0.921247805782583;-0.276851740888035;0.024331650328801;0.074968612715870;1.043269200285698;1.136150864664731];
my_HFE = -75.525761719817922;
my_energies_eV = my_energies*27.2114;
my_HFE_eV = my_HFE*27.2114;

% Their data
il_HFE = -74.961754063;
il_energies = [-20.24094;-1.27218;-0.62173;-0.45392;-0.39176;0.61293;0.75095];
il_energies_eV = il_energies*27.2114;
il_HFE_eV = il_HFE*27.2114;

% Differences
err_HFE = my_HFE - il_HFE;
err_HFE_eV = err_HFE*27.2114;
err_energies = my_energies - il_energies;
err_energies_eV = err_energies*27.2114;

% My corrected data
corrected_energies = my_energies + err_HFE;
corrected_energies_eV = corrected_energies*27.2114;
