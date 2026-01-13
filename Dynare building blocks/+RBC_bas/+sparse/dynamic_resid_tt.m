function [T_order, T] = dynamic_resid_tt(y, x, params, steady_state, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 6
    T = [T; NaN(6 - size(T, 1), 1)];
end
T(1) = y(16)^(-params(4));
T(2) = params(1)*y(30)^(-params(4));
T(3) = params(6)*exp(y(33))*(y(17)/y(32))^(params(6)-1);
T(4) = T(3)+1-params(5);
T(5) = exp(y(19))*y(3)^params(6);
T(6) = y(18)^(1-params(6));
end
