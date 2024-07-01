clearvars, clc

% Parameters
n = input("Enter matrix size: ");
v = input("Enter a value for v: ");
w = input("Enter a value for w: ");

% Matrix elements
% main diag %
M1=[0 v
    v 0];
% upper diag %
M2=[0 0
    w 0];
% lower diag %
M3=[0 w
    0 0];

finalMatrix=kron(eye(n),M1)+kron(diag(ones(n-1,1),1),M2)+kron(diag(ones(n-1,1),-1),M3);

%%% run the code below and see what matrix it can output %%%
%%% then try to understand how my code can build the same matrix like yours %%%

n=n; %%%put a number you like, better a small number so you can see how it works %%%
eye(n);
diag(ones(n-1,1),-1);
diag(ones(n-1,1),1);

eigenvalues = eig(finalMatrix);

figure(1);
plot(eigenvalues, '.')


% Sort Eigenvalues and apply indices to eigenvectors
[EV,EW] = eig(finalMatrix);
[B1, sb] = sort(real(diag(EW)));
sortedEV = EV( : , sb);

%Calculate the wave function
wavefunction = conj(sortedEV).*sortedEV;


% Find zero energies
tol = 1.e-6; % Should not check floating point numbers for exact equality, so define a tolerance
zero_energy_indices = find(abs(B1) < tol);
zero_energies = sortedEV(zero_energy_indices);
corresponding_wavefunctions = wavefunction(:, zero_energy_indices)


% Get summation of matrices
col1_corresponding_wavefunctions = corresponding_wavefunctions(:, 1)
col2_corresponding_wavefunctions = corresponding_wavefunctions(:, 2)

wavefunction_size = size(corresponding_wavefunctions(), 1);
wavefunction_half_size = wavefunction_size / 2;

sum_wavefunction1 = zeros(wavefunction_half_size, 1);
sum_wavefunction2 = zeros(wavefunction_half_size, 1);

j = 1;

for i = 1:2:wavefunction_size-1
     sum_wavefunction1(j) = col1_corresponding_wavefunctions(i) + col1_corresponding_wavefunctions(i+1);
     sum_wavefunction2(j) = col2_corresponding_wavefunctions(i) + col2_corresponding_wavefunctions(i+1);
     j = j + 1;
end

sum_wavefunction1
sum_wavefunction2


figure(2);
plot(sum_wavefunction1)
figure(3);
plot(sum_wavefunction2)
