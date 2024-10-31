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
    T = GS.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(23, 43);
g1(1,12)=1/y(19);
g1(1,35)=(-(T(2)*(T(4)/y(38)-1)));
g1(1,36)=T(2);
g1(1,19)=(-y(12))/(y(19)*y(19));
g1(1,38)=(-(T(2)*(-(y(35)*T(4)))/(y(38)*y(38))));
g1(1,39)=(-(T(2)*exp(y(39))));
g1(1,22)=(-(T(2)*(exp(y(39))-y(36)-y(35)+y(35)*T(4)/y(38))));
g1(1,24)=(-((exp(y(39))-y(36)-y(35)+y(35)*T(4)/y(38))*params(2)*exp(y(22))*(-(params(3)*exp(y(24))))));
g1(2,12)=1;
g1(2,14)=(-(getPowerDeriv(y(14),params(7),1)/T(10)));
g1(2,37)=(-((-(T(2)*getPowerDeriv(y(37),params(7),1)))/T(10)));
g1(2,22)=(-((-T(6))/T(10)));
g1(2,24)=(-((-(T(5)*params(2)*exp(y(22))*(-(params(3)*exp(y(24))))))/T(10)));
g1(3,12)=(-(params(8)*exp(y(23))*(y(17)/(1-params(3)*exp(y(24)))-1)));
g1(3,13)=1;
g1(3,17)=(-(params(8)*exp(y(23))*y(12)*1/(1-params(3)*exp(y(24)))));
g1(3,21)=(-(params(8)*exp(y(23))*exp(y(21))));
g1(3,23)=(-(T(3)+params(9)*(-(params(8)*exp(y(23))))));
g1(3,24)=(-(params(8)*exp(y(23))*y(12)*(-(y(17)*(-(params(3)*exp(y(24))))))/((1-params(3)*exp(y(24)))*(1-params(3)*exp(y(24))))));
g1(4,14)=(-1);
g1(4,1)=(-((1-params(3)*exp(y(9)))*(1-y(4))));
g1(4,15)=1;
g1(4,2)=(-((1-params(3)*exp(y(9)))*(-T(1))));
g1(4,4)=(-((1-params(3)*exp(y(9)))*(-y(1))));
g1(4,9)=(-(((1-y(4))*y(1)+T(1)*(1-y(2)))*(-(params(3)*exp(y(9))))));
g1(5,2)=(-(1-(1-params(3)*exp(y(9)))*y(3)-y(5)));
g1(5,16)=1;
g1(5,3)=(-(y(2)*(-(1-params(3)*exp(y(9))))));
g1(5,5)=(-(1-y(2)));
g1(5,9)=(-(y(2)*(-(y(3)*(-(params(3)*exp(y(9))))))));
g1(6,15)=(-(1/y(16)));
g1(6,16)=(-((-y(15))/(y(16)*y(16))));
g1(6,17)=1;
g1(7,17)=(-(T(9)*getPowerDeriv(y(17),1-params(10),1)));
g1(7,18)=1;
g1(8,17)=(-(T(9)*getPowerDeriv(y(17),(-params(10)),1)));
g1(8,19)=1;
g1(9,20)=1;
g1(9,24)=T(4)*(-(params(3)*exp(y(24))));
g1(10,6)=1;
g1(10,21)=(-1);
g1(10,25)=1;
g1(11,10)=(-1);
g1(11,25)=(-1);
g1(11,26)=1;
g1(11,11)=(-1);
g1(12,16)=(-1);
g1(12,27)=1;
g1(13,15)=(-1);
g1(13,28)=1;
g1(14,20)=(-1);
g1(14,29)=1;
g1(15,16)=(-(100*1/(steady_state(5))/(y(16)/(steady_state(5)))));
g1(15,30)=1;
g1(16,15)=(-(100*1/(steady_state(4))/(y(15)/(steady_state(4)))));
g1(16,31)=1;
g1(17,17)=(-(100*1/(steady_state(6))/(y(17)/(steady_state(6)))));
g1(17,32)=1;
g1(18,14)=(-(100*1/(steady_state(3))/(y(14)/(steady_state(3)))));
g1(18,33)=1;
g1(19,6)=(-params(11));
g1(19,21)=1;
g1(19,40)=1;
g1(20,7)=(-params(12));
g1(20,22)=1;
g1(20,41)=(-1);
g1(21,8)=(-params(13));
g1(21,23)=1;
g1(21,42)=(-1);
g1(22,9)=(-params(14));
g1(22,24)=1;
g1(22,43)=(-1);
g1(23,10)=(-1);
g1(23,34)=1;

end
