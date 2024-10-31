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
    T = flow_specification.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(32, 46);
g1(1,40)=T(1);
g1(1,13)=(-T(8))/(y(13)*y(13));
g1(1,41)=(-(T(1)*(-((1-(params(6)-params(3))/(1-params(3)))*T(8)))/(y(41)*y(41))));
g1(1,42)=(-T(1));
g1(1,22)=(-(T(11)*(1-params(3))*(-(params(2)*y(45)))/(y(22)*y(22))));
g1(1,45)=(-(T(11)*(1-params(3))*params(2)/y(22)));
g1(2,16)=(-(exp(y(28))*T(7)/T(5)));
g1(2,17)=1;
g1(2,28)=(-(T(12)/T(5)));
g1(3,8)=1;
g1(3,11)=(-(T(9)*T(8)/(1-params(3))));
g1(3,17)=(-T(9));
g1(4,9)=(-(1/y(10)));
g1(4,10)=(-((-y(9))/(y(10)*y(10))));
g1(4,11)=1;
g1(5,11)=(-(T(4)*getPowerDeriv(y(11),1-params(10),1)));
g1(5,12)=1;
g1(6,11)=(-(T(4)*getPowerDeriv(y(11),(-params(10)),1)));
g1(6,13)=1;
g1(7,10)=1;
g1(7,24)=1;
g1(8,24)=1;
g1(8,25)=(-1);
g1(8,26)=(-1);
g1(9,14)=(-(getPowerDeriv(y(14),T(2),1)));
g1(9,16)=1;
g1(10,21)=(-(getPowerDeriv(y(21),(-params(5)),1)));
g1(10,22)=1;
g1(11,18)=1;
g1(11,43)=(-(y(45)*params(2)*(1-params(3))/y(22)));
g1(11,44)=(-(y(45)*params(2)*(1-params(3))/y(22)));
g1(11,22)=(-((y(43)+y(44))*(-(y(45)*params(2)*(1-params(3))))/(y(22)*y(22))));
g1(11,45)=(-((y(43)+y(44))*params(2)*(1-params(3))/y(22)));
g1(12,16)=(-(y(26)*exp(y(28))*T(7)));
g1(12,20)=1;
g1(12,26)=(-T(12));
g1(12,28)=(-(y(26)*T(12)));
g1(13,9)=(-T(8));
g1(13,20)=1;
g1(13,21)=(-1);
g1(14,15)=1;
g1(14,25)=(-(exp(y(28))*T(7)/T(10)));
g1(14,28)=(-(exp(y(28))*y(25)*T(7)/T(10)));
g1(15,16)=(-(T(10)/T(5)));
g1(15,18)=1;
g1(16,15)=(-y(18));
g1(16,18)=(-y(15));
g1(16,20)=(-1);
g1(16,23)=1;
g1(17,14)=(-y(19));
g1(17,17)=(-y(24));
g1(17,19)=(-y(14));
g1(17,23)=1;
g1(17,24)=(-y(17));
g1(18,1)=(-(1-(1-params(3))*y(2)-params(6)));
g1(18,10)=1;
g1(18,2)=(-(y(1)*(-(1-params(3)))));
g1(19,3)=(-(1-params(3)));
g1(19,14)=1;
g1(19,4)=(-(1-params(3)));
g1(20,8)=(-(y(24)/y(23)));
g1(20,23)=(-((-(y(8)*y(24)))/(y(23)*y(23))));
g1(20,24)=(-(y(8)/y(23)));
g1(20,27)=1;
g1(21,5)=1;
g1(21,28)=(-1);
g1(21,29)=1;
g1(22,6)=(-1);
g1(22,29)=(-1);
g1(22,30)=1;
g1(22,7)=(-1);
g1(23,10)=(-1);
g1(23,31)=1;
g1(24,9)=(-1);
g1(24,32)=1;
g1(25,10)=(-(100*1/y(10)));
g1(25,33)=1;
g1(26,9)=(-(100*1/y(9)));
g1(26,34)=1;
g1(27,11)=(-(100*1/y(11)));
g1(27,35)=1;
g1(28,14)=(-(100*1/y(14)));
g1(28,36)=1;
g1(29,15)=(-(100*1/y(15)));
g1(29,37)=1;
g1(30,21)=(-(100*1/y(21)));
g1(30,38)=1;
g1(31,5)=(-params(13));
g1(31,28)=1;
g1(31,46)=1;
g1(32,6)=(-1);
g1(32,39)=1;

end
