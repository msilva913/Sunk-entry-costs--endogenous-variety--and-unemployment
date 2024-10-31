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

T(1) = params(5)/(1-params(3))/(params(6)/(1-params(3)));
T(2) = T(1)*params(4)/(params(5)+params(4));
T(3) = params(5)/(1-params(3))/T(1)^(1-params(10));
T(4) = (params(4)-params(3))/(1-params(3));
T(5) = 1-params(3)*exp(y(13));
T(6) = params(2)*exp(y(11))*T(5);
T(7) = 1-T(4);
T(8) = y(3)^params(7);
T(9) = ((T(2)-(1-params(3))*(T(2)*(1-params(6)/(1-params(3)))+T(4)*(1-params(4)/(params(5)+params(4)))))*(((1-params(2))/params(2)+params(3))/((params(1)-params(9))*params(6)*(1-params(8))/(params(6)*(1-params(8))+(1-params(2))/params(2)+params(4)+params(5)/(1-params(3))*params(8))*(1+(1-params(2))/params(2))))^(1/params(7)))^params(7);
T(10) = params(8)*exp(y(12))*(exp(y(10))-y(1)+y(1)*y(6)/T(5));

end
