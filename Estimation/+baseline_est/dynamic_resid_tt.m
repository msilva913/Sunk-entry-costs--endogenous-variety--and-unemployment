function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double  vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double  vector of endogenous variables in the order stored
%                                                    in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double  matrix of exogenous variables (in declaration order)
%                                                    for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double  vector of steady state values
%   params        [M_.param_nbr by 1]        double  vector of parameter values in declaration order
%   it_           scalar                     double  time period for exogenous variables for which
%                                                    to evaluate the model
%
% Output:
%   T           [#temp variables by 1]       double  vector of temporary terms
%

assert(length(T) >= 23);

T(1) = (1-params(2))/params(2);
T(2) = (1-params(3))*params(2)*exp(y(49))*y(77)/y(39);
T(3) = y(33)-y(22)+y(22)*y(27)/(1-params(3));
T(4) = y(74)-y(72)-y(71)+y(71)*(1-(params(6)-params(3))/(1-params(3)))/y(73);
T(5) = 1/((1+params(4))/params(4)-1);
T(6) = params(7)/(1-params(3))/(params(8)/(1-params(3)));
T(7) = params(3)*(1+T(6)*params(6)/(params(7)+params(6))-params(6)/(params(7)+params(6)));
T(8) = (1+params(4))/params(4)/((1+params(4))/params(4)-1);
T(9) = (params(3)+((1+params(4))/params(4)-1)*(T(1)+params(3)))/(params(3)+(1+params(4))/params(4)*(T(1)+params(3)));
T(10) = params(12)*y(32)*params(1)/T(8)/T(9);
T(11) = 1+1/y(29)*(T(1)+params(6))/(1-params(3));
T(12) = (y(32)*params(1)/T(8)-T(10))/T(11);
T(13) = y(32)*params(1)/T(8)-T(12)+T(6)*T(12)/(1-params(3))-params(10)*T(10);
T(14) = (T(12)*(1+T(1))/(T(1)+params(3)))^(1/params(9));
T(15) = T(7)/T(14);
T(16) = (1-params(3))*(T(8)-1)*(1-params(6)/(params(7)+params(6)))/params(13)/(T(1)+params(3)*T(8));
T(17) = exp(y(50))*(T(10)-params(10)*T(10))/T(13);
T(18) = exp(y(51))*y(24)/T(15);
T(19) = exp(y(52))*params(7)/(1-params(3))/T(6)^(1-params(11));
T(20) = y(27)^(1-params(11))*T(19);
T(21) = y(27)^(-params(11))*T(19);
T(22) = exp(y(51))*T(15)/(1+params(9));
T(23) = (y(24)/(exp(y(51))*T(15)))^(1+params(9));

end
