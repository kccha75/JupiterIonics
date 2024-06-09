% Stopping condition in ODE when concentrations become negative

function [position,isterminal,direction] = chemstopcondition(t,y)

    k=min(y);
    if k<0
        k=0;
    end
  position = k; % The value that we want to be zero
  isterminal = 1;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction
end