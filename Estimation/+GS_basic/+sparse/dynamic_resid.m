function [residual, T_order, T] = dynamic_resid(y, x, params, steady_state, T_order, T)
if nargin < 6
    T_order = -1;
    T = NaN(8, 1);
end
[T_order, T] = GS_basic.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
residual = NaN(23, 1);
    residual(1) = (y(24)/y(31)) - (T(3));
    residual(2) = (y(24)) - ((y(26)^params(7)-T(4))/T(8));
    residual(3) = (y(25)) - (T(2)+params(9)*(1-params(8)*exp(y(34))));
    residual(4) = (y(27)) - (y(26)+(1-params(3))*((1-y(8))*y(4)+T(1)*(1-y(5))));
    residual(5) = (y(28)) - (y(5)*(1-(1-params(3))*y(7))+params(4)*(1-y(5)));
    residual(6) = (y(29)) - (y(27)/y(28));
    residual(7) = (y(30)) - (y(29)^(1-params(10))*T(7));
    residual(8) = (y(31)) - (y(29)^(-params(10))*T(7));
    residual(9) = (y(35)) - (y(32)-y(9));
    residual(10) = (y(36)) - (y(35)+y(12)+y(22));
    residual(11) = (y(37)) - (log(y(25)/y(2)));
    residual(12) = (y(38)) - (y(37)+y(14)+y(23));
    residual(13) = (y(39)) - (y(28)-(steady_state(5)));
    residual(14) = (y(40)) - (y(27)-(steady_state(4)));
    residual(15) = (y(41)) - (100*log(y(28)/(steady_state(5))));
    residual(16) = (y(42)) - (100*log(y(27)/(steady_state(4))));
    residual(17) = (y(43)) - (100*log(y(29)/(steady_state(6))));
    residual(18) = (y(44)) - (100*log(y(26)/(steady_state(3))));
    residual(19) = (y(32)) - (y(9)*params(11)-x(1));
    residual(20) = (y(33)) - (params(12)*y(10)+x(2));
    residual(21) = (y(34)) - (params(13)*y(11)+x(3));
    residual(22) = (y(45)) - (y(12));
    residual(23) = (y(46)) - (y(14));
end
