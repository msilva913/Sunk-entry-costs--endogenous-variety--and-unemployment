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

assert(length(T) >= 33);

T = baseline_est.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(24) = getPowerDeriv(T(18),params(9),1);
T(25) = getPowerDeriv(y(24)/(exp(y(51))*T(15)),1+params(9),1);
T(26) = (-((y(32)*params(1)/T(8)-T(10))*(T(1)+params(6))/(1-params(3))*(-1)/(y(29)*y(29))))/(T(11)*T(11));
T(27) = exp(y(50))*(-((T(10)-params(10)*T(10))*(T(6)*T(26)/(1-params(3))-T(26))))/(T(13)*T(13));
T(28) = getPowerDeriv(T(12)*(1+T(1))/(T(1)+params(3)),1/params(9),1);
T(29) = (-(T(7)*(1+T(1))*T(26)/(T(1)+params(3))*T(28)))/(T(14)*T(14));
T(30) = params(12)*params(1)/T(8)/T(9);
T(31) = (params(1)/T(8)-T(30))/T(11);
T(32) = exp(y(50))*(T(13)*(T(30)-params(10)*T(30))-(T(10)-params(10)*T(10))*(params(1)/T(8)-T(31)+T(6)*T(31)/(1-params(3))-params(10)*T(30)))/(T(13)*T(13));
T(33) = (-(T(7)*T(28)*(1+T(1))*T(31)/(T(1)+params(3))))/(T(14)*T(14));

end
