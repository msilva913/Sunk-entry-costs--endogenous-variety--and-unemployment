function [T_order, T] = static_resid_tt(y, x, params, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 5
    T = [T; NaN(5 - size(T, 1), 1)];
end
T(1) = y(2)^(-params(4));
T(2) = params(6)*exp(y(5))*(y(3)/y(4))^(params(6)-1);
T(3) = T(2)+1-params(5);
T(4) = exp(y(5))*y(3)^params(6);
T(5) = y(4)^(1-params(6));
end
