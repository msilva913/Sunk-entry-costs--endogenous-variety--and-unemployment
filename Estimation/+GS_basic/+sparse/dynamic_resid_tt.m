function [T_order, T] = dynamic_resid_tt(y, x, params, steady_state, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 8
    T = [T; NaN(8 - size(T, 1), 1)];
end
T(1) = (params(4)-params(3))/(1-params(3));
T(2) = params(8)*exp(y(34))*(exp(y(32))-y(24)+y(24)*y(29)/(1-params(3)));
T(3) = (1-params(3))*params(2)*exp(y(33))*(exp(y(55))-y(48)-y(47)+y(47)*(1-T(1))/y(54));
T(4) = (1-params(3))*params(2)*exp(y(33))*y(49)^params(7);
T(5) = params(5)/(1-params(3))/(params(6)/(1-params(3)));
T(6) = T(5)*params(4)/(params(5)+params(4));
T(7) = params(5)/(1-params(3))/T(5)^(1-params(10));
T(8) = ((T(6)-(1-params(3))*(T(6)*(1-params(6)/(1-params(3)))+T(1)*(1-params(4)/(params(5)+params(4)))))*(((1-params(2))/params(2)+params(3))/((params(1)-params(9))*params(6)*(1-params(8))/(params(6)*(1-params(8))+(1-params(2))/params(2)+params(4)+params(5)/(1-params(3))*params(8))*(1+(1-params(2))/params(2))))^(1/params(7)))^params(7);
end
