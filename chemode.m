% Function for system of ODEs

function y=chemode(t,x)
global N2_b HA_b Li_b kap kan k1 k2 kh F R T alpha E0 k0 E j fmp fmn

% Diffusion constants
D_N2=1e1;
D_HA=1e1;
D_Li=1e1;

% Boundary layer thickness
delta_N2=0.1;
delta_HA=0.1;
delta_Li=0.1;

% % System of ODEs
y=[-kap*x(1)*x(7)+kan*x(2)+3*k2*x(4)*x(6)+2*kh*x(3)*x(6);...
    kap*x(1)*x(7)-kan*x(2)-fmp(x(2),x(3),t)+fmn(x(2),x(3),t);...
    -6*k1*x(3)*x(5)-2*kh*x(3)*x(6)+fmp(x(2),x(3),t)-fmn(x(2),x(3),t); ...
    2*k1*x(3)*x(5)-k2*x(4)*x(6);...
    -k1*x(3)*x(5)+diffusionchem(x(5),3,D_N2,delta_N2);...
    -3*k2*x(4)*x(6)-2*kh*x(3)*x(6)+diffusionchem(x(6),2,D_HA,delta_HA); ...
    -kap*x(1)*x(7)+kan*x(2)+2*kh*x(3)*x(6)+3*k2*x(4)*x(6)+diffusionchem(x(7),1,D_Li,delta_Li);...
    k2*x(4)*x(6)];

end
