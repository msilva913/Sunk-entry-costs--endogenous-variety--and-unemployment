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

assert(length(T) >= 12);

T(1) = (1-params(3))*params(2)*y(45)/y(22);
T(2) = 1/((1+params(4))/params(4)-1);
T(3) = params(7)/(1-params(3))/(params(8)/(1-params(3)));
T(4) = params(7)/(1-params(3))/T(3)^(1-params(10));
T(5) = (1+params(4))/params(4)/((1+params(4))/params(4)-1);
T(6) = params(1)/(params(11)/((params(3)+((1+params(4))/params(4)-1)*((1-params(2))/params(2)+params(3)))/(params(3)+(1+params(4))/params(4)*((1-params(2))/params(2)+params(3)))));
T(7) = T(6)*T(5)/params(12)^T(2);
T(8) = params(8)/((1-params(2))/params(2)+params(6))*(T(6)-params(1));
T(9) = (params(1)-params(9)*params(1))/(T(6)+T(3)*T(8)/(1-params(3))-params(9)*params(1));
T(10) = (1-params(3))*T(7)*(T(5)-1)*(1-params(6)/(params(7)+params(6)))/params(12)/((1-params(2))/params(2)+params(3)*T(5));
T(11) = y(42)-y(40)+(1-(params(6)-params(3))/(1-params(3)))*T(8)/y(41);
T(12) = exp(y(28))*y(16)*T(7);

end
