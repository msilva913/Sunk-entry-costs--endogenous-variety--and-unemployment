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
    T = BGM.static_resid_tt(T, y, x, params);
end
residual = zeros(19, 1);
    residual(1) = (T(1)) - (T(1)*params(1)*(1-params(3))*(y(5)+y(4))/y(5));
    residual(2) = (y(12)) - (y(3)^(1/(params(6)-1)));
    residual(3) = (y(1)) - (y(12)*y(18)*y(8));
    residual(4) = (y(9)) - (y(7)-y(8));
    residual(5) = (y(6)) - (y(12)*y(18)/(params(6)/(params(6)-1)));
    residual(6) = (y(7)) - ((T(1)*y(6)/params(4))^params(5));
    residual(7) = (y(5)) - (y(6)*params(7)/y(18));
    residual(8) = (y(4)) - (y(1)/(params(6)*y(3)));
    residual(9) = (y(1)+y(5)*y(2)) - (y(7)*y(6)+y(4)*y(3));
    residual(10) = (y(10)) - (y(1)+y(5)*y(2));
    residual(11) = (y(11)) - (y(10)/y(12));
    residual(12) = (y(3)) - ((1-params(3))*(y(3)+y(2)));
    residual(13) = (log(y(18))) - (log(y(18))*params(9)+x(1));
    residual(14) = (y(13)) - (100*log(y(10)));
    residual(15) = (y(14)) - (100*log(y(11)));
    residual(16) = (y(15)) - (100*log(y(7)));
    residual(17) = (y(16)) - (100*log(y(8)));
    residual(18) = (y(17)) - (100*log(y(9)));
    residual(19) = (y(19)) - (log(y(18))*100);

end
