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
    T = GS_basic.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(23, 43);
g1(1,13)=1/y(20);
g1(1,36)=(-(T(2)*((1-T(1))/y(39)-1)));
g1(1,37)=T(2);
g1(1,20)=(-y(13))/(y(20)*y(20));
g1(1,39)=(-(T(2)*(-(y(36)*(1-T(1))))/(y(39)*y(39))));
g1(1,40)=(-(T(2)*exp(y(40))));
g1(1,22)=(-T(9));
g1(2,13)=1;
g1(2,15)=(-(getPowerDeriv(y(15),params(7),1)/T(14)));
g1(2,38)=T(17);
g1(2,22)=(-((-T(10))/T(14)));
g1(3,13)=(-(params(8)*exp(y(23))*(y(18)/(1-params(3))-1)));
g1(3,14)=1;
g1(3,18)=(-(params(8)*exp(y(23))*y(13)*1/(1-params(3))));
g1(3,21)=(-(params(8)*exp(y(23))*exp(y(21))));
g1(3,23)=(-(T(3)+params(9)*(-(params(8)*exp(y(23))))));
g1(4,15)=(-1);
g1(4,2)=(-((1-params(3))*(1-y(5))));
g1(4,16)=1;
g1(4,3)=(-((1-params(3))*(-T(1))));
g1(4,5)=(-((1-params(3))*(-y(2))));
g1(5,3)=(-(1-(1-params(3))*y(4)-params(4)));
g1(5,17)=1;
g1(5,4)=(-(y(3)*(-(1-params(3)))));
g1(6,16)=(-(1/y(17)));
g1(6,17)=(-((-y(16))/(y(17)*y(17))));
g1(6,18)=1;
g1(7,18)=(-(T(13)*getPowerDeriv(y(18),1-params(10),1)));
g1(7,19)=1;
g1(8,18)=(-(T(13)*getPowerDeriv(y(18),(-params(10)),1)));
g1(8,20)=1;
g1(9,6)=1;
g1(9,21)=(-1);
g1(9,24)=1;
g1(10,9)=(-1);
g1(10,24)=(-1);
g1(10,25)=1;
g1(10,11)=(-1);
g1(11,1)=(-(T(15)/T(4)));
g1(11,14)=(-(T(16)/T(4)));
g1(11,26)=1;
g1(12,10)=(-1);
g1(12,26)=(-1);
g1(12,27)=1;
g1(12,12)=(-1);
g1(13,17)=(-1);
g1(13,28)=1;
g1(14,16)=(-1);
g1(14,29)=1;
g1(15,17)=(-(100*1/(steady_state(5))/T(5)));
g1(15,30)=1;
g1(16,16)=(-(100*1/(steady_state(4))/T(6)));
g1(16,31)=1;
g1(17,18)=(-(100*1/(steady_state(6))/T(7)));
g1(17,32)=1;
g1(18,15)=(-(100*1/(steady_state(3))/T(8)));
g1(18,33)=1;
g1(19,6)=(-params(11));
g1(19,21)=1;
g1(19,41)=1;
g1(20,7)=(-params(12));
g1(20,22)=1;
g1(20,42)=(-1);
g1(21,8)=(-params(13));
g1(21,23)=1;
g1(21,43)=(-1);
g1(22,9)=(-1);
g1(22,34)=1;
g1(23,10)=(-1);
g1(23,35)=1;

end
