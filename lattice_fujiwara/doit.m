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

% magnetic moment of molecule
mu = 1.5*mub;

% Lattice
lambda = 1054e-9;
kL = 2*pi/lambda;
Er = hbar^2*kL^2/(2*m);
fr= Er/h;

%% Bulk Gas Properties

fbar = 100; % trap frequency in Hz
N=1e5;
Ef = (3*N)^(1/3)*(h*fbar);
Tf = Ef/kB;
kF = sqrt(2*m*Ef/hbar^2);

%% Lattice Depth, frequency, and length scale

V0 = 200;
omega = 2*pi*sqrt(4*V0)*fr;
aho = sqrt(hbar/(mu_mass*omega));

%% Magnetic Field

% feshbach
a_bg = 166.978*a0;
Delta = 6.910;
B0 = 202.15;

% Feshbach field
B2a = @(B) a_bg*(1-Delta./(B-B0));

Bvec = linspace(100,210,1e6);
B_L = [195 207];
avec = B2a(Bvec);


%% Harmonic Contact
output = harmonic_contact;

%%
figure(10);
clf

X = aho./avec;
[X,inds] = sort(X);

% Plot the Feshbach resonance
subplot(231)
plot(Bvec,avec/a0,'linewidth',2)
ylim([-300 400])

xlim(B_L);
ylabel('scattering length $a/a_0$','interpreter','latex')
xlabel('magnetic field (G)','interpreter','latex');
set(gca,'fontsize',12,'xgrid','on','ygrid','on','fontname','times');

% Plot the feshbach resonance but in aho/a
subplot(232)
plot(Bvec,aho./avec,'linewidth',2)
% ylim([-200 200])
xlim([min(Bvec) max(Bvec)]);
ylabel('coupling $a_{\mathrm{ho}}/a$','interpreter','latex')
xlabel('magnetic field (G)','interpreter','latex');
set(gca,'fontsize',12,'xgrid','on','ygrid','on','fontname','times');
xlim(B_L);

% Energy Spectrum
subplot(234)
f1 = output(1).ainv2E;
f2 = output(2).ainv2E;
f3 = output(3).ainv2E;
plot(X,f1(X),'linewidth',2)
xlabel('$a_{\mathrm{ho}}/a$','interpreter','latex')
ylabel('energy $(\hbar\omega)$','interpreter','latex')
hold on
plot(X,f2(X),'linewidth',2)
plot(X,f3(X),'linewidth',2)
xlim([-10 10]);
set(gca,'YTick',[0 1.5 2.5 3.5 4.5 5.5],'YGrid','on','Xgrid','on')
ylim([-2 5]);
set(gca,'fontsize',12,'xgrid','on','ygrid','on','fontname','times');

% Derivative
subplot(235)
f1 = output(1).ainv2dEd1a;
f2 = output(2).ainv2dEd1a;
f3 = output(3).ainv2dEd1a;
plot(X,f1(X),'linewidth',2)
xlabel('$a_{\mathrm{ho}}/a$','interpreter','latex')
ylabel('$dE/d(1/a)~(\hbar\omega/a_{\mathrm{ho}})$','interpreter','latex')
hold on
plot(X,f2(X),'linewidth',2)
xlim([-10 10]);
ylim([-1 0]);
set(gca,'fontsize',12,'xgrid','on','ygrid','on','fontname','times');

% Contact
subplot(236)
f1 = output(1).ainv2dEd1a;
f2 = output(2).ainv2dEd1a;
f3 = output(3).ainv2dEd1a;
plot(Bvec,-8*pi*f1(aho./avec)*(a0/aho),'linewidth',2)
xlabel('magnetic field (G)','interpreter','latex');
ylabel('contact $(1/a_{0})$','interpreter','latex')
hold on
plot(Bvec,-8*pi*f2(aho./avec)*(a0/aho),'linewidth',2)
xlim([195 206]);
ylim([0 .01]);
set(gca,'fontsize',12,'xgrid','on','ygrid','on','fontname','times');

%% Output Functions

depth2aho = @(V0) sqrt(hbar/(mu_mass*(2*pi*sqrt(4*V0)*fr)));
field2contact_1 = @(B,V0) -8*pi*f1(depth2aho(V0)./B2a(B));
