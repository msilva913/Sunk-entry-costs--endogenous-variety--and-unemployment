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

assert(length(T) >= 19);

T(1) = (1-params(1))/params(1);
T(2) = 1-params(2)*exp(y(37));
T(3) = params(1)*y(58)/y(29)*T(2);
T(4) = y(17)/T(2);
T(5) = 1-(params(5)-params(2))/(1-params(2));
T(6) = 1/((1+params(3))/params(3)-1);
T(7) = params(6)/(1-params(2))/(params(7)/(1-params(2)));
T(8) = params(6)/(1-params(2))/T(7)^(1-params(11));
T(9) = (1+params(3))/params(3)/((1+params(3))/params(3)-1);
T(10) = params(14)/(params(12)/((params(2)+((1+params(3))/params(3)-1)*(T(1)+params(2)))/(params(2)+(1+params(3))/params(3)*(T(1)+params(2)))));
T(11) = T(10)*T(9)/params(13)^T(6);
T(12) = (T(10)-params(14))/(1+(T(1)+params(5))/(1-params(2))*1/(params(7)/(1-params(2))*params(8)));
T(13) = (1-params(8))/params(8)*T(12)/(params(7)/(1-params(2)));
T(14) = (params(14)-params(14)*params(10))/(T(10)-T(12)+T(7)/(1-params(2))*(T(12)+(1-params(8))/params(8)*T(12))-params(14)*params(10));
T(15) = params(2)*(1+T(7)*params(5)/(params(6)+params(5))-params(5)/(params(6)+params(5)))/(T(12)*(1+T(1))/(T(1)+params(2)))^(1/params(9));
T(16) = (1-params(2))*T(11)*(T(9)-1)*(1-params(5)/(params(6)+params(5)))/params(13)/(T(1)+params(2)*T(9));
T(17) = y(55)-y(53)-y(52)+T(5)*(y(52)/y(54)+T(13));
T(18) = exp(y(36))*y(22)*T(11);
T(19) = y(14)/T(15);

end
