function [T_order, T] = static_resid_tt(y, x, params, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 8
    T = [T; NaN(8 - size(T, 1), 1)];
end
T(1) = params(5)/(1-params(3))/(params(6)/(1-params(3)));
T(2) = T(1)*params(4)/(params(5)+params(4));
T(3) = params(5)/(1-params(3))/T(1)^(1-params(10));
T(4) = (params(4)-params(3))/(1-params(3));
T(5) = (1-params(3))*params(2)*exp(y(10))*(exp(y(9))-y(2)-y(1)+y(1)*(1-T(4))/y(8));
T(6) = y(3)^params(7);
T(7) = ((T(2)-(1-params(3))*(T(2)*(1-params(6)/(1-params(3)))+T(4)*(1-params(4)/(params(5)+params(4)))))*(((1-params(2))/params(2)+params(3))/((params(1)-params(9))*params(6)*(1-params(8))/(params(6)*(1-params(8))+(1-params(2))/params(2)+params(4)+params(5)/(1-params(3))*params(8))*(1+(1-params(2))/params(2))))^(1/params(7)))^params(7);
T(8) = params(8)*exp(y(11))*(exp(y(9))-y(1)+y(1)*y(6)/(1-params(3)));
end
