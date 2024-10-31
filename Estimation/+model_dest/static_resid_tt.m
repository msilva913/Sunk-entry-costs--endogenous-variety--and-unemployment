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

assert(length(T) >= 18);

T(1) = (1-params(1))/params(1);
T(2) = params(6)/(1-params(2))/(params(7)/(1-params(2)));
T(3) = params(6)/(1-params(2))/T(2)^(1-params(11));
T(4) = (1+params(3))/params(3)/((1+params(3))/params(3)-1);
T(5) = 1/((1+params(3))/params(3)-1);
T(6) = params(14)/(params(12)/((params(2)+((1+params(3))/params(3)-1)*(T(1)+params(2)))/(params(2)+(1+params(3))/params(3)*(T(1)+params(2)))));
T(7) = T(6)*T(4)/params(13)^T(5);
T(8) = (T(6)-params(14))/(1+(T(1)+params(5))/(1-params(2))*1/(params(7)/(1-params(2))*params(8)));
T(9) = T(8)*(1-params(8))/params(8)/(params(7)/(1-params(2)));
T(10) = (params(14)-params(14)*params(10))/(T(6)-T(8)+T(2)/(1-params(2))*(T(8)+T(8)*(1-params(8))/params(8))-params(14)*params(10));
T(11) = params(2)*(1+T(2)*params(5)/(params(6)+params(5))-params(5)/(params(6)+params(5)))/(T(8)*(1+T(1))/(T(1)+params(2)))^(1/params(9));
T(12) = (1-params(2))*T(7)*(T(4)-1)*(1-params(5)/(params(6)+params(5)))/params(13)/(T(1)+params(2)*T(4));
T(13) = T(9)+y(2)/y(9);
T(14) = 1-params(2)*exp(y(27));
T(15) = 1-(params(5)-params(2))/(1-params(2));
T(16) = exp(y(26))*T(7)*y(12);
T(17) = y(7)/T(14);
T(18) = y(4)/T(11);

end
