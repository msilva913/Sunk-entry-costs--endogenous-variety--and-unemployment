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
    T = baseline_alt.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(37, 55);
g1(1,11)=1/y(18);
g1(1,48)=(-(T(3)*((1-(params(6)-params(3))/(1-params(3)))/y(50)-1)));
g1(1,49)=T(3);
g1(1,18)=(-y(11))/(y(18)*y(18));
g1(1,50)=(-(T(3)*(-(y(48)*(1-(params(6)-params(3))/(1-params(3)))))/(y(50)*y(50))));
g1(1,51)=(-T(3));
g1(1,28)=(-(T(4)*(1-params(3))*(-(params(2)*y(54)))/(y(28)*y(28))));
g1(1,54)=(-(T(4)*(1-params(3))*params(2)/y(28)));
g1(2,16)=(-(exp(y(34))*y(21)*T(13)*(T(25)+T(26))/T(2)));
g1(2,18)=(-(T(34)/T(2)));
g1(2,21)=(-(exp(y(34))*T(14)/T(2)));
g1(2,22)=1;
g1(2,34)=(-(T(18)/T(2)));
g1(3,11)=(-(params(10)*(y(16)/(1-params(3))-1)));
g1(3,12)=1;
g1(3,16)=(-(params(10)*y(11)*1/(1-params(3))));
g1(3,22)=(-params(10));
g1(4,10)=1;
g1(4,13)=(-(1/T(17)*T(22)));
g1(4,16)=(-(T(22)*(-(y(13)*T(28)))/(T(17)*T(17))));
g1(4,18)=(-(T(22)*(-(y(13)*T(35)))/(T(17)*T(17))));
g1(5,10)=(-1);
g1(5,47)=T(3);
g1(5,11)=1;
g1(5,28)=y(47)*(1-params(3))*(-(params(2)*y(54)))/(y(28)*y(28));
g1(5,54)=y(47)*(1-params(3))*params(2)/y(28);
g1(6,14)=(-(1/y(15)));
g1(6,15)=(-((-y(14))/(y(15)*y(15))));
g1(6,16)=1;
g1(7,16)=(-(T(6)*getPowerDeriv(y(16),1-params(12),1)));
g1(7,17)=1;
g1(8,16)=(-(T(6)*getPowerDeriv(y(16),(-params(12)),1)));
g1(8,18)=1;
g1(9,15)=1;
g1(9,30)=1;
g1(10,30)=1;
g1(10,31)=(-1);
g1(10,32)=(-1);
g1(11,19)=(-(getPowerDeriv(y(19),1/(params(4)-1),1)));
g1(11,21)=1;
g1(12,26)=(-(getPowerDeriv(y(26),(-params(5)),1)));
g1(12,28)=1;
g1(13,23)=1;
g1(13,52)=(-(y(54)*params(2)*(1-params(3))/y(28)));
g1(13,53)=(-(y(54)*params(2)*(1-params(3))/y(28)));
g1(13,28)=(-((y(52)+y(53))*(-(y(54)*params(2)*(1-params(3))))/(y(28)*y(28))));
g1(13,54)=(-((y(52)+y(53))*params(2)*(1-params(3))/y(28)));
g1(14,16)=(-(y(32)*exp(y(34))*y(21)*T(13)*(T(25)+T(26))));
g1(14,18)=(-(y(32)*T(34)));
g1(14,21)=(-(y(32)*exp(y(34))*T(14)));
g1(14,25)=1;
g1(14,32)=(-T(18));
g1(14,34)=(-(y(32)*T(18)));
g1(15,25)=1;
g1(15,26)=(-1);
g1(15,27)=(-1);
g1(16,13)=(-(T(20)*1/T(17)*T(23)));
g1(16,16)=(-(T(21)*T(28)/(1+params(9))+T(20)*T(23)*(-(y(13)*T(28)))/(T(17)*T(17))));
g1(16,18)=(-(T(21)*T(35)/(1+params(9))+T(20)*T(23)*(-(y(13)*T(35)))/(T(17)*T(17))));
g1(16,27)=1;
g1(17,16)=(-((T(15)*exp(y(34))*y(31)*T(13)*(T(25)+T(26))-exp(y(34))*y(31)*T(14)*T(29))/(T(15)*T(15))));
g1(17,18)=(-((T(15)*exp(y(34))*y(31)*T(33)-exp(y(34))*y(31)*T(14)*T(36))/(T(15)*T(15))));
g1(17,20)=1;
g1(17,31)=(-(exp(y(34))*T(14)/T(15)));
g1(17,34)=(-(exp(y(34))*y(31)*T(14)/T(15)));
g1(18,16)=(-(y(21)*T(29)/T(2)));
g1(18,18)=(-(y(21)*T(36)/T(2)));
g1(18,21)=(-(T(15)/T(2)));
g1(18,23)=1;
g1(19,20)=(-y(23));
g1(19,23)=(-y(20));
g1(19,25)=(-1);
g1(19,29)=1;
g1(20,19)=(-y(24));
g1(20,22)=(-y(30));
g1(20,24)=(-y(19));
g1(20,29)=1;
g1(20,30)=(-y(22));
g1(21,13)=(-1);
g1(21,1)=(-((1-params(3))*(1-y(4))));
g1(21,14)=1;
g1(21,2)=(-((1-params(3))*(-((params(6)-params(3))/(1-params(3))))));
g1(21,4)=(-((1-params(3))*(-y(1))));
g1(22,2)=(-(1-(1-params(3))*y(3)-params(6)));
g1(22,15)=1;
g1(22,3)=(-(y(2)*(-(1-params(3)))));
g1(23,5)=(-(1-params(3)));
g1(23,19)=1;
g1(23,6)=(-(1-params(3)));
g1(24,12)=(-(y(30)/y(29)));
g1(24,29)=(-((-(y(12)*y(30)))/(y(29)*y(29))));
g1(24,30)=(-(y(12)/y(29)));
g1(24,33)=1;
g1(25,7)=1;
g1(25,34)=(-1);
g1(25,35)=1;
g1(26,8)=(-1);
g1(26,35)=(-1);
g1(26,36)=1;
g1(26,9)=(-1);
g1(27,15)=(-1);
g1(27,37)=1;
g1(28,14)=(-1);
g1(28,38)=1;
g1(29,15)=(-(100*1/y(15)));
g1(29,39)=1;
g1(30,14)=(-(100*1/y(14)));
g1(30,40)=1;
g1(31,16)=(-(100*1/y(16)));
g1(31,41)=1;
g1(32,13)=(-(100*1/y(13)));
g1(32,42)=1;
g1(33,19)=(-(100*1/y(19)));
g1(33,43)=1;
g1(34,20)=(-(100*1/y(20)));
g1(34,44)=1;
g1(35,26)=(-(100*1/y(26)));
g1(35,45)=1;
g1(36,7)=(-params(14));
g1(36,34)=1;
g1(36,55)=1;
g1(37,8)=(-1);
g1(37,46)=1;

end
