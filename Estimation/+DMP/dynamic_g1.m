function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g1
%

if T_flag
    T = DMP.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(10, 17);
g1(1,14)=params(2);
g1(1,9)=(-T(5))/(y(9)*y(9));
g1(1,15)=(-(params(2)*(-((1-params(4))*T(5)))/(y(15)*y(15))));
g1(1,16)=(-(params(2)*exp(y(16))));
g1(2,4)=1;
g1(2,7)=(-(params(7)*T(5)));
g1(2,10)=(-(params(7)*exp(y(10))));
g1(3,1)=(-(1-y(2)-params(4)));
g1(3,6)=1;
g1(3,2)=y(1);
g1(4,5)=(-(1/y(6)));
g1(4,6)=(-((-y(5))/(y(6)*y(6))));
g1(4,7)=1;
g1(5,7)=(-(T(4)*getPowerDeriv(y(7),1-params(9),1)));
g1(5,8)=1;
g1(6,7)=(-(T(4)*getPowerDeriv(y(7),(-params(9)),1)));
g1(6,9)=1;
g1(7,6)=(-(100*1/(steady_state(3))/T(1)));
g1(7,11)=1;
g1(8,5)=(-(100*1/(steady_state(2))/T(2)));
g1(8,12)=1;
g1(9,7)=(-(100*1/(steady_state(4))/T(3)));
g1(9,13)=1;
g1(10,3)=(-params(10));
g1(10,10)=1;
g1(10,17)=1;

end
