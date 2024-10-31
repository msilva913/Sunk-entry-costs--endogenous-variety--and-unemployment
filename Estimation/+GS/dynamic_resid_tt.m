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

assert(length(T) >= 10);

T(1) = (params(4)-params(3))/(1-params(3));
T(2) = params(2)*exp(y(22))*(1-params(3)*exp(y(24)));
T(3) = params(8)*exp(y(23))*(exp(y(21))-y(12)+y(12)*y(17)/(1-params(3)*exp(y(24))));
T(4) = 1-T(1);
T(5) = y(37)^params(7);
T(6) = T(2)*T(5);
T(7) = params(5)/(1-params(3))/(params(6)/(1-params(3)));
T(8) = T(7)*params(4)/(params(5)+params(4));
T(9) = params(5)/(1-params(3))/T(7)^(1-params(10));
T(10) = ((T(8)-(1-params(3))*(T(8)*(1-params(6)/(1-params(3)))+T(1)*(1-params(4)/(params(5)+params(4)))))*(((1-params(2))/params(2)+params(3))/((params(1)-params(9))*params(6)*(1-params(8))/(params(6)*(1-params(8))+(1-params(2))/params(2)+params(4)+params(5)/(1-params(3))*params(8))*(1+(1-params(2))/params(2))))^(1/params(7)))^params(7);

end
