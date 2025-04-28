function [residual, T_order, T] = static_resid(y, x, params, T_order, T)
if nargin < 5
    T_order = -1;
    T = NaN(8, 1);
end
[T_order, T] = GS_basic.sparse.static_resid_tt(y, x, params, T_order, T);
residual = NaN(23, 1);
    residual(1) = (y(1)/y(8)) - (T(5));
    residual(2) = (y(1)) - ((T(6)-(1-params(3))*params(2)*exp(y(10))*T(6))/T(7));
    residual(3) = (y(2)) - (T(8)+params(9)*(1-params(8)*exp(y(11))));
    residual(4) = (y(4)) - (y(3)+(1-params(3))*(y(4)*(1-y(8))+T(4)*(1-y(5))));
    residual(5) = (y(5)) - (y(5)*(1-(1-params(3))*y(7))+params(4)*(1-y(5)));
    residual(6) = (y(6)) - (y(4)/y(5));
    residual(7) = (y(7)) - (T(3)*y(6)^(1-params(10)));
    residual(8) = (y(8)) - (T(3)*y(6)^(-params(10)));
residual(9) = y(12);
    residual(10) = (y(13)) - (y(12)+y(12)+y(22));
residual(11) = y(14);
    residual(12) = (y(15)) - (y(14)+y(14)+y(23));
    residual(13) = (y(16)) - (y(5)-(y(5)));
    residual(14) = (y(17)) - (y(4)-(y(4)));
    residual(15) = (y(18)) - (100*log(y(5)/(y(5))));
    residual(16) = (y(19)) - (100*log(y(4)/(y(4))));
    residual(17) = (y(20)) - (100*log(y(6)/(y(6))));
    residual(18) = (y(21)) - (100*log(y(3)/(y(3))));
    residual(19) = (y(9)) - (y(9)*params(11)-x(1));
    residual(20) = (y(10)) - (y(10)*params(12)+x(2));
    residual(21) = (y(11)) - (y(11)*params(13)+x(3));
    residual(22) = (y(22)) - (y(12));
    residual(23) = (y(23)) - (y(14));
end
