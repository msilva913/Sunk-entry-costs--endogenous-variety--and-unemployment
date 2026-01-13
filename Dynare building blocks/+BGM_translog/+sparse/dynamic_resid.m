function [residual, T_order, T] = dynamic_resid(y, x, params, steady_state, T_order, T)
if nargin < 6
    T_order = -1;
    T = NaN(3, 1);
end
[T_order, T] = BGM_translog.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
residual = NaN(21, 1);
    residual(1) = (y(41)) - (1+1/(params(9)*y(24)));
    residual(2) = (y(42)) - (y(41)/(y(41)-1));
    residual(3) = (y(33)) - (sqrt(y(24))*T(1));
    residual(4) = (T(2)) - (T(3)*(y(47)+y(46))/y(26));
    residual(5) = (y(22)) - (y(33)*y(39)*y(29));
    residual(6) = (y(27)) - (y(33)*y(39)/y(41));
    residual(7) = (y(28)) - ((T(2)*y(27)/params(4))^params(5));
    residual(8) = (y(30)) - (y(28)-y(29));
    residual(9) = (y(26)) - (y(27)*params(6)/y(39));
    residual(10) = (y(25)) - ((y(41)-1)/y(41)*y(22)/y(24));
    residual(11) = (y(22)+y(26)*y(23)) - (y(27)*y(28)+y(24)*y(25));
    residual(12) = (y(31)) - (y(22)+y(26)*y(23));
    residual(13) = (y(32)) - (y(31)/y(33));
    residual(14) = (y(24)) - ((1-params(3))*(y(3)+y(2)));
    residual(15) = (log(y(39))) - (params(8)*log(y(18))+x(1));
    residual(16) = (y(34)) - (100*log(y(31)));
    residual(17) = (y(35)) - (100*log(y(32)));
    residual(18) = (y(36)) - (100*log(y(28)));
    residual(19) = (y(37)) - (100*log(y(29)));
    residual(20) = (y(38)) - (100*log(y(30)));
    residual(21) = (y(40)) - (log(y(39))*100);
end
