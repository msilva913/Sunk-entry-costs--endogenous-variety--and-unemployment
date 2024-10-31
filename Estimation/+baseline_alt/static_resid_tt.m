function T = static_resid_tt(T, y, x, params)
% function T = static_resid_tt(T, y, x, params)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%
% Output:
%   T         [#temp variables by 1]  double   vector of temporary terms
%

assert(length(T) >= 19);

T(1) = (1-params(2))/params(2);
T(2) = params(4)/(params(4)-1);
T(3) = params(7)/(1-params(3))/(params(8)/(1-params(3)));
T(4) = params(7)/(1-params(3))/T(3)^(1-params(12));
T(5) = params(3)*(1+T(3)*params(6)/(params(7)+params(6))-params(6)/(params(7)+params(6)));
T(6) = 1/y(9)*(T(1)+params(6))/(1-params(3));
T(7) = (params(10)*y(7)+(T(1)+params(6))/y(9))/((1-params(3))*(1-params(10)));
T(8) = T(6)/T(7);
T(9) = (params(1)-params(11)*params(1)*T(8))/(1-T(8))-params(11)*params(1);
T(10) = T(9)/T(7);
T(11) = T(2)/params(13)^(1/(params(4)-1));
T(12) = ((params(1)-params(11)*params(1)*T(8))/(1-T(8))+T(10))*T(11);
T(13) = (1-params(3))*(T(2)-1)*(1-params(6)/(params(7)+params(6)))*T(12)/params(13)/(T(1)+T(2)*params(3));
T(14) = (T(10)*(1+T(1))/(T(1)+params(3)))^(1/params(9));
T(15) = T(5)/T(14);
T(16) = exp(y(25))*T(12)*y(12);
T(17) = y(4)/T(15);
T(18) = T(15)/(1+params(9));
T(19) = T(17)^(1+params(9));

end
