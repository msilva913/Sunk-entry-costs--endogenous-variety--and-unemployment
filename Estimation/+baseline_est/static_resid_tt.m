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

assert(length(T) >= 22);

T(1) = (1-params(2))/params(2);
T(2) = params(7)/(1-params(3))/(params(8)/(1-params(3)));
T(3) = params(3)*(1+T(2)*params(6)/(params(7)+params(6))-params(6)/(params(7)+params(6)));
T(4) = (1+params(4))/params(4)/((1+params(4))/params(4)-1);
T(5) = 1/((1+params(4))/params(4)-1);
T(6) = (params(3)+((1+params(4))/params(4)-1)*(T(1)+params(3)))/(params(3)+(1+params(4))/params(4)*(T(1)+params(3)));
T(7) = y(12)*params(1)/T(4)*params(12)/T(6);
T(8) = 1+1/y(9)*(T(1)+params(6))/(1-params(3));
T(9) = (y(12)*params(1)/T(4)-T(7))/T(8);
T(10) = y(12)*params(1)/T(4)-T(9)+T(2)*T(9)/(1-params(3))-T(7)*params(10);
T(11) = (T(9)*(1+T(1))/(T(1)+params(3)))^(1/params(9));
T(12) = T(3)/T(11);
T(13) = (1-params(3))*(T(4)-1)*(1-params(6)/(params(7)+params(6)))/params(13)/(T(1)+params(3)*T(4));
T(14) = (1-params(3))*params(2)*exp(y(29))*(y(13)-y(3)-y(2)+y(2)*(1-(params(6)-params(3))/(1-params(3)))/y(9));
T(15) = y(13)-y(2)+y(2)*y(7)/(1-params(3));
T(16) = (T(7)-T(7)*params(10))/T(10)*exp(y(30));
T(17) = exp(y(31))*y(4)/T(12);
T(18) = params(7)/(1-params(3))/T(2)^(1-params(11))*exp(y(32));
T(19) = y(7)^(1-params(11))*T(18);
T(20) = T(18)*y(7)^(-params(11));
T(21) = T(12)*exp(y(31))/(1+params(9));
T(22) = (y(4)/(T(12)*exp(y(31))))^(1+params(9));

end
