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

assert(length(T) >= 10);

T(1) = params(7)/(1-params(3))/(params(8)/(1-params(3)));
T(2) = params(7)/(1-params(3))/T(1)^(1-params(10));
T(3) = (1+params(4))/params(4)/((1+params(4))/params(4)-1);
T(4) = 1/((1+params(4))/params(4)-1);
T(5) = params(1)/(params(11)/((params(3)+((1+params(4))/params(4)-1)*((1-params(2))/params(2)+params(3)))/(params(3)+(1+params(4))/params(4)*((1-params(2))/params(2)+params(3)))));
T(6) = T(5)*T(3)/params(12)^T(4);
T(7) = params(8)/((1-params(2))/params(2)+params(6))*(T(5)-params(1));
T(8) = (params(1)-params(9)*params(1))/(T(5)+T(1)*T(7)/(1-params(3))-params(9)*params(1));
T(9) = (1-params(3))*T(6)*(T(3)-1)*(1-params(6)/(params(7)+params(6)))/params(12)/((1-params(2))/params(2)+params(3)*T(3));
T(10) = exp(y(21))*T(6)*y(9);

end
