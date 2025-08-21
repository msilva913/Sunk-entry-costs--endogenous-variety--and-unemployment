function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = model_dest.static_resid_tt(T, y, x, params);
end
residual = zeros(41, 1);
    residual(1) = (T(13)) - (params(1)*T(14)*(y(13)-y(3)-y(2)+T(13)*T(15)));
    residual(2) = (y(13)) - (T(16)/T(4));
    residual(3) = (y(3)) - (T(10)*(y(13)-y(2)+T(17)*(y(2)+T(9)*y(9)))+params(14)*params(10)*(1-T(10)));
    residual(4) = (y(1)) - (T(18)^params(9));
    residual(5) = (y(2)) - (y(1)-params(1)*T(14)*y(1));
    residual(6) = (y(7)) - (y(5)/y(6));
    residual(7) = (y(8)) - (T(3)*y(7)^(1-params(11)));
    residual(8) = (y(9)) - (T(3)*y(7)^(-params(11)));
    residual(9) = (y(21)) - (1-y(6));
    residual(10) = (y(21)) - (y(23)+y(22));
    residual(11) = (y(12)) - (y(10)^T(5));
    residual(12) = (y(19)) - (y(17)^(-params(4)));
    residual(13) = (y(14)) - (params(1)*T(14)*(y(14)+y(15)));
    residual(14) = (y(16)) - (T(16)*y(23));
    residual(15) = (y(16)) - (y(17)+y(18)+y(9)*T(9)*y(5));
    residual(16) = (y(18)) - (T(11)/(1+params(9))*T(18)^(1+params(9)));
    residual(17) = (y(11)) - (exp(y(26))*T(7)*y(22)/T(12));
    residual(18) = (y(14)) - (T(12)*y(12)/T(4));
    residual(19) = (y(20)) - (y(16)+y(14)*y(11));
    residual(20) = (y(20)) - (y(13)*y(21)+y(10)*y(15));
    residual(21) = (y(5)) - (y(4)+T(14)*(y(5)*(1-y(9))+(params(5)-params(2))/(1-params(2))*(1-y(6))));
    residual(22) = (y(6)) - (y(6)*(1-T(14)*y(8))+(1-y(6))*y(25));
    residual(23) = (y(10)) - (T(14)*(y(10)+y(11)));
    residual(24) = (y(24)) - (y(3)*y(21)/y(20));
    residual(25) = (y(25)) - (1-T(14)*T(15));
residual(26) = y(28);
    residual(27) = (y(29)) - (y(28)+y(28)+y(41));
    residual(28) = (y(30)) - (y(6)-(y(6)));
    residual(29) = (y(31)) - (y(5)-(y(5)));
    residual(30) = (y(32)) - (100*log(y(6)));
    residual(31) = (y(33)) - (100*log(y(5)));
    residual(32) = (y(34)) - (100*log(y(7)));
    residual(33) = (y(35)) - (100*log(y(4)));
    residual(34) = (y(36)) - (100*log(y(10)));
    residual(35) = (y(37)) - (100*log(y(11)));
    residual(36) = (y(39)) - (100*log(y(17)));
    residual(37) = (y(40)) - (100*log(y(20)));
    residual(38) = (y(38)) - (100*log(y(12)));
    residual(39) = (y(26)) - (y(26)*params(15)-x(1));
    residual(40) = (y(27)) - (y(27)*params(16)+x(2));
    residual(41) = (y(41)) - (y(28));

end
