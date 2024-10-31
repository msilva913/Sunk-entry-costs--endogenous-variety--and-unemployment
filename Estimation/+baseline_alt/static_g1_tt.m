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

assert(length(T) >= 34);

T = baseline_alt.static_resid_tt(T, y, x, params);

T(20) = getPowerDeriv(T(17),params(9),1);
T(21) = getPowerDeriv(T(17),1+params(9),1);
T(22) = (-(T(6)*params(10)/((1-params(3))*(1-params(10)))))/(T(7)*T(7));
T(23) = ((1-T(8))*(-(params(11)*params(1)*T(22)))-(params(1)-params(11)*params(1)*T(8))*(-T(22)))/((1-T(8))*(1-T(8)));
T(24) = (T(7)*T(23)-T(9)*params(10)/((1-params(3))*(1-params(10))))/(T(7)*T(7));
T(25) = getPowerDeriv(T(10)*(1+T(1))/(T(1)+params(3)),1/params(9),1);
T(26) = (-(T(5)*(1+T(1))*T(24)/(T(1)+params(3))*T(25)))/(T(14)*T(14));
T(27) = (1-params(3))*(T(2)-1)*(1-params(6)/(params(7)+params(6)))*T(11)*(T(23)+T(24))/params(13)/(T(1)+T(2)*params(3));
T(28) = (-(T(1)+params(6)))/(y(9)*y(9))/((1-params(3))*(1-params(10)));
T(29) = (T(7)*(T(1)+params(6))/(1-params(3))*(-1)/(y(9)*y(9))-T(6)*T(28))/(T(7)*T(7));
T(30) = ((1-T(8))*(-(params(11)*params(1)*T(29)))-(params(1)-params(11)*params(1)*T(8))*(-T(29)))/((1-T(8))*(1-T(8)));
T(31) = T(11)*(T(30)+(T(7)*T(30)-T(9)*T(28))/(T(7)*T(7)));
T(32) = exp(y(25))*y(12)*T(31);
T(33) = (-(T(5)*T(25)*(1+T(1))*(T(7)*T(30)-T(9)*T(28))/(T(7)*T(7))/(T(1)+params(3))))/(T(14)*T(14));
T(34) = (1-params(3))*(T(2)-1)*(1-params(6)/(params(7)+params(6)))*T(31)/params(13)/(T(1)+T(2)*params(3));

end
