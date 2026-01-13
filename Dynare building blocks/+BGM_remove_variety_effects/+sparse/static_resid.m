function [residual, T_order, T] = static_resid(y, x, params, T_order, T)
if nargin < 5
    T_order = -1;
    T = NaN(1, 1);
end
[T_order, T] = BGM_remove_variety_effects.sparse.static_resid_tt(y, x, params, T_order, T);
residual = NaN(19, 1);
    residual(1) = (T(1)) - (T(1)*params(1)*(1-params(3))*(y(5)+y(4))/y(5));
    residual(2) = (y(12)) - (1.0);
    residual(3) = (y(1)) - (y(8)*1.0*y(18));
    residual(4) = (y(9)) - (y(7)-y(8));
    residual(5) = (y(6)) - (1.0*y(18)/(params(6)/(params(6)-1)));
    residual(6) = (y(7)) - ((T(1)*y(6)/params(4))^params(5));
    residual(7) = (y(5)) - (y(6)*params(7)/y(18));
    residual(8) = (y(4)) - (y(1)/(params(6)*y(3)));
    residual(9) = (y(1)+y(5)*y(2)) - (y(7)*y(6)+y(4)*y(3));
    residual(10) = (y(10)) - (y(1)+y(5)*y(2));
    residual(11) = (y(11)) - (y(10)/1.0);
    residual(12) = (y(3)) - ((1-params(3))*(y(3)+y(2)));
    residual(13) = (log(y(18))) - (log(y(18))*params(9)+x(1));
    residual(14) = (y(13)) - (100*log(y(10)));
    residual(15) = (y(14)) - (100*log(y(11)));
    residual(16) = (y(15)) - (100*log(y(7)));
    residual(17) = (y(16)) - (100*log(y(8)));
    residual(18) = (y(17)) - (100*log(y(9)));
    residual(19) = (y(19)) - (log(y(18))*100);
end
