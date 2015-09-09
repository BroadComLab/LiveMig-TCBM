function [value]=delta_f(arg)
    % The function return the delta dirac
    if arg==0 
        value=1;
    else
        value=0;
    end
end