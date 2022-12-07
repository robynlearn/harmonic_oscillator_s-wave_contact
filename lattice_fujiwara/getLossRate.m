function [R1, R2] = getLossRate(Vin,Bin,Omega,gamma)
% Fundamental Constants
amu = 1.66054e-27; % amu in kg
h = 6.62607015e-34; % planck's constant  in Js
hbar = h/(2*pi); % reduced planck's constant Js
a0 = 5.29177e-11;       % bohr radisu in m
mubh = 1.39962449e6; % bohr magneton/h in Hz/Gauss
mub = mubh*h; % (Bohr magneton in J/Gauss).
kB = 1.381e-23 ; % boltzmann constant in J/K

% magnetic moment of molecule
mu = 1.5*mub;

% atom mass
m = 40*amu; % amss
mu_mass = m/2;

% Lattice
lambda = 1054e-9;
kL = 2*pi/lambda;
Er = hbar^2*kL^2/(2*m);
fr= Er/h;

%% Contact
[C1,C2] = getHarmonicContact(Vin,Bin);

%% d(1/a)/dB
% feshbach
a_bg = 166.978*a0;
Delta = 6.910;
B0 = 202.15;

% Feshbach field
B2a = @(B) a_bg*(1-Delta./(B-B0));

dainvdB = (1-a_bg./B2a(Bin)).^2./(a_bg*Delta);

%% Rate

R1 = (Omega^2/gamma)*(C1*hbar^2/(4*pi*m)*(1/mu)*dainvdB);
R2 = (Omega^2/gamma)*(C2*hbar^2/(4*pi*m)*(1/mu)*dainvdB);

end

