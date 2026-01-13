function [residual, T_order, T] = static_resid(y, x, params, T_order, T)
if nargin < 5
    T_order = -1;
    T = NaN(2, 1);
end
[T_order, T] = BGM_translog.sparse.static_resid_tt(y, x, params, T_order, T);
residual = NaN(21, 1);
    residual(1) = (y(20)) - (1+1/(params(9)*y(3)));
    residual(2) = (y(21)) - (y(20)/(y(20)-1));
    residual(3) = (y(12)) - (sqrt(y(3))*T(1));
    residual(4) = (T(2)) - (T(2)*params(1)*(1-params(3))*(y(5)+y(4))/y(5));
    residual(5) = (y(1)) - (y(12)*y(18)*y(8));
    residual(6) = (y(6)) - (y(12)*y(18)/y(20));
    residual(7) = (y(7)) - ((T(2)*y(6)/params(4))^params(5));
    residual(8) = (y(9)) - (y(7)-y(8));
    residual(9) = (y(5)) - (y(6)*params(6)/y(18));
    residual(10) = (y(4)) - ((y(20)-1)/y(20)*y(1)/y(3));
    residual(11) = (y(1)+y(5)*y(2)) - (y(6)*y(7)+y(3)*y(4));
    residual(12) = (y(10)) - (y(1)+y(5)*y(2));
    residual(13) = (y(11)) - (y(10)/y(12));
    residual(14) = (y(3)) - ((1-params(3))*(y(3)+y(2)));
    residual(15) = (log(y(18))) - (log(y(18))*params(8)+x(1));
    residual(16) = (y(13)) - (100*log(y(10)));
    residual(17) = (y(14)) - (100*log(y(11)));
    residual(18) = (y(15)) - (100*log(y(7)));
    residual(19) = (y(16)) - (100*log(y(8)));
    residual(20) = (y(17)) - (100*log(y(9)));
    residual(21) = (y(19)) - (log(y(18))*100);
end
