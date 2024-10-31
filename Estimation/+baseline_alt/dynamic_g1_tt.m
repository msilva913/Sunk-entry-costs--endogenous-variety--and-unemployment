function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
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

assert(length(T) >= 36);

T = baseline_alt.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(22) = getPowerDeriv(T(19),params(9),1);
T(23) = getPowerDeriv(T(19),1+params(9),1);
T(24) = (-(T(8)*params(10)/((1-params(3))*(1-params(10)))))/(T(9)*T(9));
T(25) = ((1-T(10))*(-(params(11)*params(1)*T(24)))-(params(1)-params(11)*params(1)*T(10))*(-T(24)))/((1-T(10))*(1-T(10)));
T(26) = (T(9)*T(25)-T(11)*params(10)/((1-params(3))*(1-params(10))))/(T(9)*T(9));
T(27) = getPowerDeriv(T(12)*(1+T(1))/(T(1)+params(3)),1/params(9),1);
T(28) = (-(T(7)*(1+T(1))*T(26)/(T(1)+params(3))*T(27)))/(T(16)*T(16));
T(29) = (1-params(3))*(T(2)-1)*(1-params(6)/(params(7)+params(6)))*T(13)*(T(25)+T(26))/params(13)/(T(1)+T(2)*params(3));
T(30) = (-(T(1)+params(6)))/(y(18)*y(18))/((1-params(3))*(1-params(10)));
T(31) = (T(9)*(T(1)+params(6))/(1-params(3))*(-1)/(y(18)*y(18))-T(8)*T(30))/(T(9)*T(9));
T(32) = ((1-T(10))*(-(params(11)*params(1)*T(31)))-(params(1)-params(11)*params(1)*T(10))*(-T(31)))/((1-T(10))*(1-T(10)));
T(33) = T(13)*(T(32)+(T(9)*T(32)-T(11)*T(30))/(T(9)*T(9)));
T(34) = exp(y(34))*y(21)*T(33);
T(35) = (-(T(7)*T(27)*(1+T(1))*(T(9)*T(32)-T(11)*T(30))/(T(9)*T(9))/(T(1)+params(3))))/(T(16)*T(16));
T(36) = (1-params(3))*(T(2)-1)*(1-params(6)/(params(7)+params(6)))*T(33)/params(13)/(T(1)+T(2)*params(3));

end
