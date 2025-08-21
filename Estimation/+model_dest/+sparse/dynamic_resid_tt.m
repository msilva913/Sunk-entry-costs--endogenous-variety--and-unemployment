function [T_order, T] = dynamic_resid_tt(y, x, params, steady_state, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 19
    T = [T; NaN(19 - size(T, 1), 1)];
end
T(1) = (1-params(1))/params(1);
T(2) = 1-params(2)*exp(y(68));
T(3) = params(1)*y(101)/y(60)*T(2);
T(4) = y(48)/T(2);
T(5) = 1-(params(5)-params(2))/(1-params(2));
T(6) = 1/((1+params(3))/params(3)-1);
T(7) = params(6)/(1-params(2))/(params(7)/(1-params(2)));
T(8) = params(6)/(1-params(2))/T(7)^(1-params(11));
T(9) = (1+params(3))/params(3)/((1+params(3))/params(3)-1);
T(10) = params(14)/(params(12)/((params(2)+((1+params(3))/params(3)-1)*(T(1)+params(2)))/(params(2)+(1+params(3))/params(3)*(T(1)+params(2)))));
T(11) = T(10)*T(9)/params(13)^T(6);
T(12) = (T(10)-params(14))/(1+(T(1)+params(5))/(1-params(2))*1/(params(7)/(1-params(2))*params(8)));
T(13) = (1-params(8))/params(8)*T(12)/(params(7)/(1-params(2)));
T(14) = (params(14)-params(14)*params(10))/(T(10)-T(12)+T(7)/(1-params(2))*(T(12)+(1-params(8))/params(8)*T(12))-params(14)*params(10));
T(15) = params(2)*(1+T(7)*params(5)/(params(6)+params(5))-params(5)/(params(6)+params(5)))/(T(12)*(1+T(1))/(T(1)+params(2)))^(1/params(9));
T(16) = (1-params(2))*T(11)*(T(9)-1)*(1-params(5)/(params(6)+params(5)))/params(13)/(T(1)+params(2)*T(9));
T(17) = y(95)-y(85)-y(84)+T(5)*(y(84)/y(91)+T(13));
T(18) = exp(y(67))*y(53)*T(11);
T(19) = y(45)/T(15);
end
