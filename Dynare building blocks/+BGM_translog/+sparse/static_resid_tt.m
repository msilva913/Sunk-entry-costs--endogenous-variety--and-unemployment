function [T_order, T] = static_resid_tt(y, x, params, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 2
    T = [T; NaN(2 - size(T, 1), 1)];
end
T(1) = exp(1/(y(3)*2*params(9)));
T(2) = y(1)^(-params(7));
end
