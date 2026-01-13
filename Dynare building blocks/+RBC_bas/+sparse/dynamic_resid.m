function [residual, T_order, T] = dynamic_resid(y, x, params, steady_state, T_order, T)
if nargin < 6
    T_order = -1;
    T = NaN(6, 1);
end
[T_order, T] = RBC_bas.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
residual = NaN(14, 1);
    residual(1) = (T(1)) - (T(2)*T(4));
    residual(2) = (params(2)*y(18)^(1/params(3))) - (T(1)*y(21));
    residual(3) = (y(17)) - ((1-params(5))*y(3)+y(22));
    residual(4) = (y(15)) - (y(16)+y(22));
    residual(5) = (y(15)) - (T(5)*T(6));
    residual(6) = (y(21)) - (y(15)*(1-params(6))/y(18));
    residual(7) = (y(20)) - (y(15)*params(6)*4/y(3));
    residual(8) = (y(19)) - (params(7)*y(5)+x(1));
    residual(9) = (y(23)) - (100*log(y(15)));
    residual(10) = (y(24)) - (100*log(y(17)));
    residual(11) = (y(25)) - (100*log(y(16)));
    residual(12) = (y(26)) - (100*log(y(18)));
    residual(13) = (y(27)) - (100*log(y(21)));
    residual(14) = (y(28)) - (100*log(y(22)));
end
