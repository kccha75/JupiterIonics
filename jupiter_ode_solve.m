clear;clc;%close all
global N2_b HA_b Li_b kap kan k1 k2 kh F R T alpha E0 k0 E j fmp fmn

% Solves system of ODEs reflecting chemical equations

% Final time
tf=100;

% ODE=0 (small time) or algebraic equation=1 (large time)
ode=0;

% Constants for BV equation
F=9.65*10^4;    % Faraday's constant
R=8.31 ;        % Gas constant
T=298.2;        % Temperature

alpha=1/2;
E0=3.04;        % equilibrium Voltage

k0=10^-4;   

% BV equation functions
E=@(t) abs(t-50)/50+3;
j=@(Li_p,Li_0,t) k0*(Li_p.*exp(alpha*F/(R*T).*(E(t)-E0)))-k0*(Li_0.*exp(-(1-alpha)*F/(R*T)*(E(t)-E0)));

fmp=@(Li_p,Li_0,t) k0*(Li_p.*exp(alpha*F/(R*T).*(E(t)-E0)));
fmn=@(Li_p,Li_0,t) k0*(Li_0.*exp(-(1-alpha)*F/(R*T)*(E(t)-E0)));

% Reference in ODE solver
%
% S             x(1)
% Li+ads        x(2)
% Li0ads        x(3)
% Li3N          x(4)
% N2            x(5)
% HA            x(6)
% Li+           x(7)
% NH3           x(8)

% Initial concentrations at x=0 t=0
S_0=10;
Liads_0=0;
Li0ads_0=0;
Li3N_0=0;
N2_0=10;
HA_0=10;
Li_0=100;
NH3_0=0;

% Initial bulk concentrations
N2_b=10;
HA_b=10;
Li_b=100;

% Reaction rates
kap=2;
kan=1;
k1=4;
k2=10;
kh=10;

% ODE setup for algebraic solver
M=eye(8);
if ode==0
    M(5,:)=0;
    M(6,:)=0;
    M(7,:)=0;
end

% Set ODE tolerance, stopping condition and algebraic equations
odeoptions=odeset('RelTol', 1e-7, 'AbsTol',1e-7,'Events',@chemstopcondition,'Mass',M);

% Stiff ODE solver 
[t,y]=ode15s(@chemode,[0,tf],[S_0,Liads_0,Li0ads_0,Li3N_0,N2_0,HA_0,Li_0,NH3_0],odeoptions);

% Current calculation
jj=j(y(:,2),y(:,3),t);

% Plots vs time
figure('Position',[200 300 600 800]); fsz=15; lw=2;

subplot(5,2,1)
plot(t,E(t),'linewidth',lw)
ylabel('$E$','interpreter','latex','fontsize',fsz)

subplot(5,2,2)
plot(t,jj,'linewidth',lw)
ylabel('$j$','interpreter','latex','fontsize',fsz)

subplot(5,2,3)
plot(t,y(:,1),'linewidth',lw)
ylabel('$S^*$','interpreter','latex','fontsize',fsz)

subplot(5,2,4)
plot(t,y(:,2),'linewidth',lw)
ylabel('$\mathrm{Li}^+_{(\mathrm{ads})}$','interpreter','latex','fontsize',fsz)

subplot(5,2,5)
plot(t,y(:,3),'linewidth',lw)
ylabel('$\mathrm{Li}^0_{(\mathrm{ads})}$','interpreter','latex','fontsize',fsz)

subplot(5,2,6)
plot(t,y(:,4),'linewidth',lw)
ylabel('$\mathrm{Li}_3\mathrm{N}$','interpreter','latex','fontsize',fsz)

subplot(5,2,7)
plot(t,y(:,5),'linewidth',lw)
ylabel('$\mathrm{N}_2$','interpreter','latex','fontsize',fsz)

subplot(5,2,8)
plot(t,y(:,6),'linewidth',lw)
ylabel('$\mathrm{HA}$','interpreter','latex','fontsize',fsz)

subplot(5,2,9)
plot(t,y(:,7),'linewidth',lw)
ylabel('$\mathrm{Li}^+_{(\mathrm{sol})}$','interpreter','latex','fontsize',fsz)
xlabel('$t$','interpreter','latex','fontsize',fsz)

subplot(5,2,10)
plot(t,y(:,8),'linewidth',lw)
ylabel('$\mathrm{NH}_3$','interpreter','latex','fontsize',fsz)
xlabel('$t$','interpreter','latex','fontsize',fsz)

