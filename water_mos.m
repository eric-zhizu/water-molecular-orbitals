%% Import data from online text files
clear
% Data can be found on http://vergil.chemistry.gatech.edu/h2oints.txt
% 1) Geometric coordinates of the 3 nuclei in H20
% 2) Overlap integrals
% 3) 1-electron integrals
% 4) 2-electron integrals
fileID = fopen('water_minimal_basis_data.txt');
geometry = cell2mat(textscan(fileID, '%*s %f %f %f', 'HeaderLines',23));
geometry = geometry/0.529177249; % convert from Angstroms to Bohrs 
fileID = fopen('water_minimal_basis_data.txt');
overlap = cell2mat(textscan(fileID, '%*d %f %f %f %f %f %f %f', 'HeaderLines',40));
fileID = fopen('water_minimal_basis_data.txt');
one_electron = cell2mat(textscan(fileID, '%*d %f %f %f %f %f %f %f', 'HeaderLines',82));
fileID = fopen('water_minimal_basis_data.txt');
two_electron = cell2mat(textscan(fileID, '%*s %f %f %*s %f %f) %*s %f', 'HeaderLines',95));
two_electron(:,1:4) = two_electron(:,1:4) + 1; % to correct AO number

% Data can be found on http://chemistry.illinoisstate.edu/standard/che460/handouts/460water.pdf
% 5) Initial guess coefficients for LCAO
fileID = fopen('water_matrix_data.txt');
c = cell2mat(textscan(fileID, '%*d %*s %*s %f %f %f %f %f', 'HeaderLines',6));

fclose('all');

% NOTE: The potential energy matrix and overlap matrix in the Georgia Tech document is
% slightly different than the potential energy matrix and overlap matrix in the Illinois State document.
% As a result, the one_electron matrix = kinetic + potential is slightly
% different in both documents, which may have resulted in a small percent
% error in my results

%% Beginning the SCF iterations

HFE = zeros(1,31);
HFE(1) = 0; % setting Hartree-Fock Energy to 0 at the beginning for arbitrary comparison
for k=1:30 % max number of iterations
%% Density matrix P

c = c(:,1:5); % to cut off unoccupied orbitals

P = zeros(7,7);
for i=1:7
    for j=1:7
        P(i,j) = 2*sum(conj(c(i,:)).*c(j,:));
    end
end

%% Fock matrix F
F = zeros(7,7); % Fock matrix
f = zeros(7,7); % Two electron integral component of Fock matrix
for i=1:7
    for j=1:7
        for l=1:7
            for s=1:7
                % The two electron integral data was given in a certain
                % order where for any [pq|rs] p > q > r > s
                % The following code sorts the variables provided in the
                % conditions of the for loop. I use this to search the
                % imported data two_electron for the two-electron integral
                % values
                te_one = zeros(1,4); 
                if i >= j
                    te_one(1) = i;
                    te_one(2) = j;
                else
                    te_one(1) = j;
                    te_one(2) = i;
                end
                if l >= s
                    te_one(3) = l;
                    te_one(4) = s;
                else
                    te_one(3) = s;
                    te_one(4) = l;
                end
                if te_one(3) > te_one(1)
                    te_one = te_one([3 4 1 2]);
                elseif te_one(3) == te_one(1)
                    if te_one(4) > te_one(2)
                        te_one = te_one([3 4 1 2]);
                    end
                end
                coulomb = 0;
                for row = 1:length(two_electron)
                    if(two_electron(row,1:4) == te_one)
                        coulomb = two_electron(row,5);
                        break
                    end
                end

                te_two = zeros(1,4);
                if i >= l
                    te_two(1) = i;
                    te_two(2) = l;
                else
                    te_two(1) = l;
                    te_two(2) = i;
                end
                if j >= s
                    te_two(3) = j;
                    te_two(4) = s;
                else
                    te_two(3) = s;
                    te_two(4) = j;
                end
                if te_two(3) > te_two(1)
                    te_two = te_two([3 4 1 2]);
                elseif te_two(3) == te_two(1)
                    if te_two(4) > te_two(2)
                        te_two = te_two([3 4 1 2]);
                    end
                end
                exchange = 0;
                for row = 1:length(two_electron)
                    if(two_electron(row,1:4) == te_two)
                        exchange = two_electron(row,5);
                        break
                    end
                end
                
                % I do a double summation of the two-electron integrals
                % by constructing a matrix size l x s, and summing up the
                % entire matrix f(:)
                f(l,s) = P(l,s)*(coulomb - .5*exchange);
            end
        end
        F(i,j) = one_electron(i,j) + sum(f(:));
    end
end

%% Calculate the energies of the orbitals
% because the secular determinant = 0

syms E
energies = double(vpasolve(det(F - E*overlap)==0,E));

%% Solve for coefficients
% Due to non-zero overlap, my matrix equation is not an eigenvector problem
% I use singular value decomposition to solve (F-energies*overlap)(c) = 0
c = zeros(7,7);
for i=1:7
    [U,S,V] = svd(F-energies(i)*overlap);
    c(:,i) = V(:,end);
end

%% Hartree-Fock Energy
expectation = sum(energies + diag(one_electron)); % <E>
repulsion_ab = 16/sqrt(sum((geometry(1,:)-geometry(2,:)).^2));
repulsion_bc = 1/sqrt(sum((geometry(2,:)-geometry(3,:)).^2));
repulsion_ca = 16/sqrt(sum((geometry(1,:)-geometry(2,:)).^2));
repulsion = repulsion_ab + repulsion_bc + repulsion_ca; % nuclear repulsion energy

HFE(k+1) = expectation + repulsion;

% I check to see if the change in energy between the two iterations is less
% than 1e-8 hartrees
if abs((HFE(k+1)-HFE(k))) < 1e-8
    break
end
end

disp(k) % number of iterations
format long
disp(HFE(HFE~=0)') % Hartree-Fock energies of each iteration
format short
disp(HFE(k+1)-HFE(k)) % difference in the last two Hartree-Fock energies
disp(energies); % MO energies
disp(c); % AO coefficients

%% STO-3G
% I grab my atomic orbitals from the STO-3G basis set. Since I now have the
% coefficients to construct my LCAO MOs, I can plot my orbitals as pretty
% pictures
fileID = fopen('water_sto3g.txt');
h_1s = cell2mat(textscan(fileID, '%f %f', 'HeaderLines',16));
fileID = fopen('water_sto3g.txt');
o_1s = cell2mat(textscan(fileID, '%f %f', 'HeaderLines',22));
fileID = fopen('water_sto3g.txt');
o_2sp = cell2mat(textscan(fileID, '%f %f %f', 'HeaderLines',26));

fclose('all');

%% Coordinate System
% I use the code from pset5, orbitals_part3_redacted as a model to graph my
% molecular orbitals for water

% First, create the coordinate axes
[Y, X, Z] = meshgrid(linspace(-3.2,3.2,100), linspace(-3.2,3.2,100), linspace(-3.2,3.2,100));

R_O = (X.^2 + Y.^2 + (Z-geometry(1,3)).^2).^.5;
R_H1 = (X.^2 + (Y-geometry(2,2)).^2 + (Z-geometry(2,3)).^2).^.5;
R_H2 = (X.^2 + (Y-geometry(3,2)).^2 + (Z-geometry(3,3)).^2).^.5;

Rperp_O = (X.^2 + Y.^2).^.5;
Rperp_H1 = (X.^2 + (Y-geometry(2,2)).^2).^.5;
Rperp_H2 = (X.^2 + (Y-geometry(3,2)).^2).^.5;

theta_O = atan2(Rperp_O, Z-geometry(1,3)); 
theta_H1 = atan2(Rperp_H1, Z-geometry(2,3)); 
theta_H2 = atan2(Rperp_H2, Z-geometry(3,3));

phi_O = atan2(Y, X);
phi_H1 = atan2(Y-geometry(2,2), X);
phi_H2 = atan2(Y-geometry(3,2), X);


%% STO-3G
% Second, create the atomic orbitals centered around each nucleus
% p. 181 of Szabo and Ostlund
H1_1s = zeros(size(R_H1));
for i=1:3
    H1_1s = H1_1s + h_1s(i,2)*((8*h_1s(i,1)^3/pi^3)^(1/4)*exp(-h_1s(i,1)*R_H1.^2));
end

H2_1s = zeros(size(R_H2));
for i=1:3
    H2_1s = H2_1s + h_1s(i,2)*((8*h_1s(i,1)^3/pi^3)^(1/4)*exp(-h_1s(i,1)*R_H2.^2));
end

O_1s = zeros(size(R_O));
for i=1:3
    O_1s = O_1s + o_1s(i,2)*((8*o_1s(i,1)^3/pi^3)^(1/4)*exp(-o_1s(i,1)*R_O.^2));
end

O_2s = zeros(size(R_O));
for i=1:3
    O_2s = O_2s + o_2sp(i,2)*((8*o_2sp(i,1)^3/pi^3)^(1/4)*exp(-o_2sp(i,1)*R_O.^2));
end

O_2px = zeros(size(R_O));
for i=1:3
    O_2px = O_2px + o_2sp(i,3)*((128*o_2sp(i,1)^5/pi^3)^(1/4)*R_O.*sin(theta_O).*cos(phi_O).*exp(-o_2sp(i,1)*R_O.^2));
end

O_2py = zeros(size(R_O));
for i=1:3
    O_2py = O_2py + o_2sp(i,3)*((128*o_2sp(i,1)^5/pi^3)^(1/4)*R_O.*sin(theta_O).*sin(phi_O).*exp(-o_2sp(i,1)*R_O.^2));
end

O_2pz = zeros(size(R_O));
for i=1:3
    O_2pz = O_2pz + o_2sp(i,3)*((128*o_2sp(i,1)^5/pi^3)^(1/4)*R_O.*cos(theta_O).*exp(-o_2sp(i,1)*R_O.^2));
end

% Normalize each orbital
H1_1s = H1_1s/sqrt(sum(H1_1s(:).^2));
H2_1s = H2_1s/sqrt(sum(H2_1s(:).^2));
O_1s = O_1s/sqrt(sum(O_1s(:).^2));
O_2s = O_2s/sqrt(sum(O_2s(:).^2));
O_2px = O_2px/sqrt(sum(O_2px(:).^2));
O_2py = O_2py/sqrt(sum(O_2py(:).^2));
O_2pz = O_2pz/sqrt(sum(O_2pz(:).^2));

% In retrospect, I realize that I could have used these orbitals to
% generate my one_electron hamiltonian matrix and my overlap matrix
% instead of having to import them from online. For example, I could have
% done the following to produce my one-electron integrals (this is pseudocode):
% for i = 1:7
%     for j = 1:7
%         one_electron(i,j) = conj(AO_i(:))'*[transform matrix]*AO_j(:);
%     end
% end
% However, constructing those transform matrices (summation of laplacians
% and potential energy matrices) was something that was beyond my abilities
% Even though the overlap matrix doesn't use any transform matrix, I felt
% that it was easier to be consistent - instead of using an overlap matrix
% I generated, I just imported the data from online as well.

% Construct molecular orbitals
MOs = matlab.lang.makeUniqueStrings(repmat({'MO'}, 1, 7), 'MO');
for i=1:7
    MOs{i} = c(1,i)*O_1s+c(2,i)*O_2s+c(3,i)*O_2px+c(4,i)*O_2py+c(5,i)*O_2pz+c(6,i)*H1_1s+c(7,i)*H2_1s;
    MOs{i} = MOs{i}/sqrt(sum(MOs{i}(:).^2)); % normalize
end

%% Artificial Correction of Error

standard_MO5 = -.39176; % energy of non-bonding orbital calculated in Dr. Jean Standard's document
correction = standard_MO5 - energies(5);
correction_eV = correction*27.2114;
corr_E = energies + correction;
corr_E_eV = corr_E*27.2114; % In my MO graphics, I use the corrected energies

%% Generate Pretty Pictures
% plotiso is a function created by Adam Cohen (2011) for pset5

figure(1)
for i=1:7
    subplot(2,4,i)
    plotiso(X, Y, Z, MOs{i}, -.1*max(MOs{i}(:)), 'green');
    plotiso(X, Y, Z, MOs{i}, .1*max(MOs{i}(:)), 'yellow');
    if i <= 5
        occ = 'Occupied ';
    else
        occ = 'Unoccupied ';
    end
    title([occ 'Orbital ' num2str(i) ': Energy = ' num2str(corr_E_eV(i)) ' eV'])
    plotiso(X,Y,Z,R_H1,.2,'cyan');
    plotiso(X,Y,Z,R_H2,.2,'cyan');
    plotiso(X,Y,Z,R_O,.2,'red');
    hold on
    plot3(geometry(1:2,1),geometry(1:2,2),geometry(1:2,3),'k:','LineWidth',1);
    plot3(geometry([1 3],1),geometry([1 3],2),geometry([1 3],3),'k:','LineWidth',1);
    if i == 5
        view([90 90])
    else
        view([90 0])
    end
end
hold off