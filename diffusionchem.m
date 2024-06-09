% Diffusion estimated to be D/dx(C_bulk-C_x=0) for some chemical C
%
% C - parameter to determine which compound C involved
% D - diffusion constant
% C_bulk - bulk concentration for compound C
% C_x - surface concentration for compount C

function D3=diffusionchem(C_x,C,D,dx)
global Li_b HA_b N2_b

if C==1
    C_b=Li_b;
elseif C==2
    C_b=HA_b;
elseif C==3
    C_b=N2_b;
end

D3=D/dx*(C_b-C_x);

end