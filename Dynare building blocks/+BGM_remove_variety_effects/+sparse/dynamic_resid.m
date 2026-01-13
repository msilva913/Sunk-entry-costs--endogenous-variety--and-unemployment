function [residual, T_order, T] = dynamic_resid(y, x, params, steady_state, T_order, T)
if nargin < 6
    T_order = -1;
    T = NaN(2, 1);
end
[T_order, T] = BGM_remove_variety_effects.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
residual = NaN(19, 1);
    residual(1) = (T(1)) - (T(2)*(y(43)+y(42))/y(24));
    residual(2) = (y(31)) - (1.0);
    residual(3) = (y(20)) - (y(27)*1.0*y(37));
    residual(4) = (y(28)) - (y(26)-y(27));
    residual(5) = (y(25)) - (1.0*y(37)/(params(6)/(params(6)-1)));
    residual(6) = (y(26)) - ((T(1)*y(25)/params(4))^params(5));
    residual(7) = (y(24)) - (y(25)*params(7)/y(37));
    residual(8) = (y(23)) - (y(20)/(params(6)*y(22)));
    residual(9) = (y(20)+y(24)*y(21)) - (y(26)*y(25)+y(23)*y(22));
    residual(10) = (y(29)) - (y(20)+y(24)*y(21));
    residual(11) = (y(30)) - (y(29)/1.0);
    residual(12) = (y(22)) - ((1-params(3))*(y(3)+y(2)));
    residual(13) = (log(y(37))) - (params(9)*log(y(18))+x(1));
    residual(14) = (y(32)) - (100*log(y(29)));
    residual(15) = (y(33)) - (100*log(y(30)));
    residual(16) = (y(34)) - (100*log(y(26)));
    residual(17) = (y(35)) - (100*log(y(27)));
    residual(18) = (y(36)) - (100*log(y(28)));
    residual(19) = (y(38)) - (log(y(37))*100);
end
