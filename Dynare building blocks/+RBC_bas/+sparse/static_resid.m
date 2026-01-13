function [residual, T_order, T] = static_resid(y, x, params, T_order, T)
if nargin < 5
    T_order = -1;
    T = NaN(5, 1);
end
[T_order, T] = RBC_bas.sparse.static_resid_tt(y, x, params, T_order, T);
residual = NaN(14, 1);
    residual(1) = (T(1)) - (T(1)*params(1)*T(3));
    residual(2) = (params(2)*y(4)^(1/params(3))) - (T(1)*y(7));
    residual(3) = (y(3)) - (y(3)*(1-params(5))+y(8));
    residual(4) = (y(1)) - (y(2)+y(8));
    residual(5) = (y(1)) - (T(4)*T(5));
    residual(6) = (y(7)) - (y(1)*(1-params(6))/y(4));
    residual(7) = (y(6)) - (y(1)*params(6)*4/y(3));
    residual(8) = (y(5)) - (y(5)*params(7)+x(1));
    residual(9) = (y(9)) - (100*log(y(1)));
    residual(10) = (y(10)) - (100*log(y(3)));
    residual(11) = (y(11)) - (100*log(y(2)));
    residual(12) = (y(12)) - (100*log(y(4)));
    residual(13) = (y(13)) - (100*log(y(7)));
    residual(14) = (y(14)) - (100*log(y(8)));
end
