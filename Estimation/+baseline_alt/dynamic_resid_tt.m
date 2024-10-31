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

assert(length(T) >= 21);

T(1) = (1-params(2))/params(2);
T(2) = params(4)/(params(4)-1);
T(3) = (1-params(3))*params(2)*y(54)/y(28);
T(4) = y(51)-y(49)-y(48)+y(48)*(1-(params(6)-params(3))/(1-params(3)))/y(50);
T(5) = params(7)/(1-params(3))/(params(8)/(1-params(3)));
T(6) = params(7)/(1-params(3))/T(5)^(1-params(12));
T(7) = params(3)*(1+T(5)*params(6)/(params(7)+params(6))-params(6)/(params(7)+params(6)));
T(8) = 1/y(18)*(T(1)+params(6))/(1-params(3));
T(9) = (params(10)*y(16)+(T(1)+params(6))/y(18))/((1-params(3))*(1-params(10)));
T(10) = T(8)/T(9);
T(11) = (params(1)-params(11)*params(1)*T(10))/(1-T(10))-params(11)*params(1);
T(12) = T(11)/T(9);
T(13) = T(2)/params(13)^(1/(params(4)-1));
T(14) = ((params(1)-params(11)*params(1)*T(10))/(1-T(10))+T(12))*T(13);
T(15) = (1-params(3))*(T(2)-1)*(1-params(6)/(params(7)+params(6)))*T(14)/params(13)/(T(1)+T(2)*params(3));
T(16) = (T(12)*(1+T(1))/(T(1)+params(3)))^(1/params(9));
T(17) = T(7)/T(16);
T(18) = exp(y(34))*y(21)*T(14);
T(19) = y(13)/T(17);
T(20) = T(17)/(1+params(9));
T(21) = T(19)^(1+params(9));

end
