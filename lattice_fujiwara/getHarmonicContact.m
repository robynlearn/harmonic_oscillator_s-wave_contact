function [C1,C2] = getHarmonicContact(Vin,Bin)
%% Constants
% Fundamental Constants
amu = 1.66054e-27; % amu in kg
h = 6.62607015e-34; % planck's constant  in Js
hbar = h/(2*pi); % reduced planck's constant Js
a0 = 5.29177e-11;       % bohr radisu in m
mubh = 1.39962449e6; % bohr magneton/h in Hz/Gauss
mub = mubh*h; % (Bohr magneton in J/Gauss).
kB = 1.381e-23 ; % boltzmann constant in J/K

% atom mass
m = 40*amu; % amss
mu_mass = m/2;

% Lattice
lambda = 1054e-9;
kL = 2*pi/lambda;
Er = hbar^2*kL^2/(2*m);
fr= Er/h;

%% Harmonic Approximation
omega = 2*pi*sqrt(4*Vin)*fr;

% Harmonic oscillator lengthscale
aho = sqrt(hbar/(mu_mass*omega));

%% Magnetic Field

% feshbach
a_bg = 166.978*a0;
Delta = 6.910;
B0 = 202.15;

% Feshbach field
B2a = @(B) a_bg*(1-Delta./(B-B0));

% Scattering Length
a = B2a(Bin);

%% Normalize scattering length to harmonic length

X = aho./a;

%% Harmonic Contact
output = harmonic_contact;

%% Get contact

f1 = output(1).ainv2dEd1a;
f2 = output(2).ainv2dEd1a;

C1 = -8*pi*f1(X)*(a0./aho);
C2 = -8*pi*f2(X)*(a0./aho);

end

