function [value]=sum_f(vect, from, to, step)
    % The function return the delta dirac
    value = 0;
    for i=from:to
        value = value + vect(i*step+1+1);
    end
end