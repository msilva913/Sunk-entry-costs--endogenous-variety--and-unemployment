function T = static_g1_tt(T, y, x, params)
% function T = static_g1_tt(T, y, x, params)
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

assert(length(T) >= 32);

T = baseline_est.static_resid_tt(T, y, x, params);

T(23) = getPowerDeriv(T(17),params(9),1);
T(24) = getPowerDeriv(y(4)/(T(12)*exp(y(31))),1+params(9),1);
T(25) = (-((y(12)*params(1)/T(4)-T(7))*(T(1)+params(6))/(1-params(3))*(-1)/(y(9)*y(9))))/(T(8)*T(8));
T(26) = exp(y(30))*(-((T(7)-T(7)*params(10))*(T(2)*T(25)/(1-params(3))-T(25))))/(T(10)*T(10));
T(27) = getPowerDeriv(T(9)*(1+T(1))/(T(1)+params(3)),1/params(9),1);
T(28) = (-(T(3)*(1+T(1))*T(25)/(T(1)+params(3))*T(27)))/(T(11)*T(11));
T(29) = params(12)*params(1)/T(4)/T(6);
T(30) = (params(1)/T(4)-T(29))/T(8);
T(31) = exp(y(30))*(T(10)*(T(29)-params(10)*T(29))-(T(7)-T(7)*params(10))*(params(1)/T(4)-T(30)+T(2)*T(30)/(1-params(3))-params(10)*T(29)))/(T(10)*T(10));
T(32) = (-(T(3)*T(27)*(1+T(1))*T(30)/(T(1)+params(3))))/(T(11)*T(11));

end
