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
    T = model_dest.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(41, 61);
g1(1,12)=1/y(19);
g1(1,53)=(-(T(3)*(T(5)*1/y(55)-1)));
g1(1,54)=T(3);
g1(1,19)=(-y(12))/(y(19)*y(19));
g1(1,55)=(-(T(3)*T(5)*(-y(53))/(y(55)*y(55))));
g1(1,56)=(-T(3));
g1(1,29)=(-(T(17)*T(2)*(-(params(1)*y(59)))/(y(29)*y(29))));
g1(1,59)=(-(T(17)*T(2)*params(1)/y(29)));
g1(1,37)=(-(T(17)*params(1)*y(59)/y(29)*(-(params(2)*exp(y(37))))));
g1(2,22)=(-(exp(y(36))*T(11)/T(9)));
g1(2,23)=1;
g1(2,36)=(-(T(18)/T(9)));
g1(3,12)=(-(T(14)*(T(4)-1)));
g1(3,13)=1;
g1(3,17)=(-(T(14)*(y(12)+y(19)*T(13))*1/T(2)));
g1(3,19)=(-(T(14)*T(4)*T(13)));
g1(3,23)=(-T(14));
g1(3,37)=(-(T(14)*(y(12)+y(19)*T(13))*(-(y(17)*(-(params(2)*exp(y(37))))))/(T(2)*T(2))));
g1(4,11)=1;
g1(4,14)=(-(1/T(15)*getPowerDeriv(T(19),params(9),1)));
g1(5,11)=(-1);
g1(5,52)=T(3);
g1(5,12)=1;
g1(5,29)=y(52)*T(2)*(-(params(1)*y(59)))/(y(29)*y(29));
g1(5,59)=y(52)*T(2)*params(1)/y(29);
g1(5,37)=y(52)*params(1)*y(59)/y(29)*(-(params(2)*exp(y(37))));
g1(6,15)=(-(1/y(16)));
g1(6,16)=(-((-y(15))/(y(16)*y(16))));
g1(6,17)=1;
g1(7,17)=(-(T(8)*getPowerDeriv(y(17),1-params(11),1)));
g1(7,18)=1;
g1(8,17)=(-(T(8)*getPowerDeriv(y(17),(-params(11)),1)));
g1(8,19)=1;
g1(9,16)=1;
g1(9,31)=1;
g1(10,31)=1;
g1(10,32)=(-1);
g1(10,33)=(-1);
g1(11,20)=(-(getPowerDeriv(y(20),T(6),1)));
g1(11,22)=1;
g1(12,27)=(-(getPowerDeriv(y(27),(-params(4)),1)));
g1(12,29)=1;
g1(13,24)=1;
g1(13,57)=(-(y(59)*params(1)*T(2)/y(29)));
g1(13,58)=(-(y(59)*params(1)*T(2)/y(29)));
g1(13,29)=(-((y(57)+y(58))*(-(y(59)*params(1)*T(2)))/(y(29)*y(29))));
g1(13,59)=(-((y(57)+y(58))*params(1)*T(2)/y(29)));
g1(13,37)=(-((y(57)+y(58))*y(59)*params(1)*(-(params(2)*exp(y(37))))/y(29)));
g1(14,22)=(-(y(33)*exp(y(36))*T(11)));
g1(14,26)=1;
g1(14,33)=(-T(18));
g1(14,36)=(-(y(33)*T(18)));
g1(15,15)=(-(y(19)*T(13)));
g1(15,19)=(-(y(15)*T(13)));
g1(15,26)=1;
g1(15,27)=(-1);
g1(15,28)=(-1);
g1(16,14)=(-(T(15)/(1+params(9))*1/T(15)*getPowerDeriv(T(19),1+params(9),1)));
g1(16,28)=1;
g1(17,21)=1;
g1(17,32)=(-(exp(y(36))*T(11)/T(16)));
g1(17,36)=(-(exp(y(36))*y(32)*T(11)/T(16)));
g1(18,22)=(-(T(16)/T(9)));
g1(18,24)=1;
g1(19,21)=(-y(24));
g1(19,24)=(-y(21));
g1(19,26)=(-1);
g1(19,30)=1;
g1(20,20)=(-y(25));
g1(20,23)=(-y(31));
g1(20,25)=(-y(20));
g1(20,30)=1;
g1(20,31)=(-y(23));
g1(21,14)=(-1);
g1(21,1)=(-(T(2)*(1-y(4))));
g1(21,15)=1;
g1(21,2)=(-(T(2)*(-((params(5)-params(2))/(1-params(2))))));
g1(21,4)=(-(T(2)*(-y(1))));
g1(21,37)=(-(((1-y(4))*y(1)+(params(5)-params(2))/(1-params(2))*(1-y(2)))*(-(params(2)*exp(y(37))))));
g1(22,2)=(-(1-T(2)*y(3)-y(35)));
g1(22,16)=1;
g1(22,3)=(-(y(2)*(-T(2))));
g1(22,35)=(-(1-y(2)));
g1(22,37)=(-(y(2)*(-(y(3)*(-(params(2)*exp(y(37))))))));
g1(23,5)=(-T(2));
g1(23,20)=1;
g1(23,6)=(-T(2));
g1(23,37)=(-((y(5)+y(6))*(-(params(2)*exp(y(37))))));
g1(24,13)=(-(y(31)/y(30)));
g1(24,30)=(-((-(y(13)*y(31)))/(y(30)*y(30))));
g1(24,31)=(-(y(13)/y(30)));
g1(24,34)=1;
g1(25,35)=1;
g1(25,37)=T(5)*(-(params(2)*exp(y(37))));
g1(26,7)=1;
g1(26,36)=(-1);
g1(26,38)=1;
g1(27,9)=(-1);
g1(27,38)=(-1);
g1(27,39)=1;
g1(27,10)=(-1);
g1(28,16)=(-1);
g1(28,40)=1;
g1(29,15)=(-1);
g1(29,41)=1;
g1(30,16)=(-(100*1/y(16)));
g1(30,42)=1;
g1(31,15)=(-(100*1/y(15)));
g1(31,43)=1;
g1(32,17)=(-(100*1/y(17)));
g1(32,44)=1;
g1(33,14)=(-(100*1/y(14)));
g1(33,45)=1;
g1(34,20)=(-(100*1/y(20)));
g1(34,46)=1;
g1(35,21)=(-(100*1/y(21)));
g1(35,47)=1;
g1(36,27)=(-(100*1/y(27)));
g1(36,49)=1;
g1(37,30)=(-(100*1/y(30)));
g1(37,50)=1;
g1(38,22)=(-(100*1/y(22)));
g1(38,48)=1;
g1(39,7)=(-params(15));
g1(39,36)=1;
g1(39,60)=1;
g1(40,8)=(-params(16));
g1(40,37)=1;
g1(40,61)=(-1);
g1(41,9)=(-1);
g1(41,51)=1;

end
