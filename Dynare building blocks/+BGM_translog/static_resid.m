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
    T = BGM_translog.static_resid_tt(T, y, x, params);
end
residual = zeros(21, 1);
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
