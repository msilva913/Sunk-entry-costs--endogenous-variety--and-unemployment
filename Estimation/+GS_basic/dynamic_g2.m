function g2 = dynamic_g2(T, y, x, params, steady_state, it_, T_flag)
% function g2 = dynamic_g2(T, y, x, params, steady_state, it_, T_flag)
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
%   g2
%

if T_flag
    T = GS_basic.dynamic_g2_tt(T, y, x, params, steady_state, it_);
end
g2_i = zeros(48,1);
g2_j = zeros(48,1);
g2_v = zeros(48,1);

g2_i(1)=1;
g2_i(2)=1;
g2_i(3)=1;
g2_i(4)=1;
g2_i(5)=1;
g2_i(6)=1;
g2_i(7)=1;
g2_i(8)=1;
g2_i(9)=1;
g2_i(10)=1;
g2_i(11)=1;
g2_i(12)=1;
g2_i(13)=1;
g2_i(14)=1;
g2_i(15)=1;
g2_i(16)=1;
g2_i(17)=2;
g2_i(18)=2;
g2_i(19)=2;
g2_i(20)=2;
g2_i(21)=2;
g2_i(22)=3;
g2_i(23)=3;
g2_i(24)=3;
g2_i(25)=3;
g2_i(26)=3;
g2_i(27)=3;
g2_i(28)=3;
g2_i(29)=3;
g2_i(30)=3;
g2_i(31)=3;
g2_i(32)=4;
g2_i(33)=4;
g2_i(34)=5;
g2_i(35)=5;
g2_i(36)=6;
g2_i(37)=6;
g2_i(38)=6;
g2_i(39)=7;
g2_i(40)=8;
g2_i(41)=11;
g2_i(42)=11;
g2_i(43)=11;
g2_i(44)=11;
g2_i(45)=15;
g2_i(46)=16;
g2_i(47)=17;
g2_i(48)=18;
g2_j(1)=536;
g2_j(2)=830;
g2_j(3)=1544;
g2_j(4)=1670;
g2_j(5)=1527;
g2_j(6)=939;
g2_j(7)=1570;
g2_j(8)=940;
g2_j(9)=837;
g2_j(10)=1673;
g2_j(11)=1656;
g2_j(12)=942;
g2_j(13)=1717;
g2_j(14)=1699;
g2_j(15)=943;
g2_j(16)=925;
g2_j(17)=617;
g2_j(18)=1629;
g2_j(19)=1613;
g2_j(20)=941;
g2_j(21)=925;
g2_j(22)=534;
g2_j(23)=744;
g2_j(24)=539;
g2_j(25)=959;
g2_j(26)=754;
g2_j(27)=964;
g2_j(28)=881;
g2_j(29)=883;
g2_j(30)=967;
g2_j(31)=969;
g2_j(32)=48;
g2_j(33)=174;
g2_j(34)=90;
g2_j(35)=132;
g2_j(36)=662;
g2_j(37)=704;
g2_j(38)=705;
g2_j(39)=749;
g2_j(40)=749;
g2_j(41)=1;
g2_j(42)=14;
g2_j(43)=560;
g2_j(44)=573;
g2_j(45)=705;
g2_j(46)=661;
g2_j(47)=749;
g2_j(48)=617;
g2_v(1)=(-1)/(y(20)*y(20));
g2_v(2)=g2_v(1);
g2_v(3)=(-(T(2)*(-(1-T(1)))/(y(39)*y(39))));
g2_v(4)=g2_v(3);
g2_v(5)=(-(T(2)*((1-T(1))/y(39)-1)));
g2_v(6)=g2_v(5);
g2_v(7)=T(2);
g2_v(8)=g2_v(7);
g2_v(9)=(-((-y(13))*(y(20)+y(20))))/(y(20)*y(20)*y(20)*y(20));
g2_v(10)=(-(T(2)*(-((-(y(36)*(1-T(1))))*(y(39)+y(39))))/(y(39)*y(39)*y(39)*y(39))));
g2_v(11)=(-(T(2)*(-(y(36)*(1-T(1))))/(y(39)*y(39))));
g2_v(12)=g2_v(11);
g2_v(13)=(-(T(2)*exp(y(40))));
g2_v(14)=(-(T(2)*exp(y(40))));
g2_v(15)=g2_v(14);
g2_v(16)=(-T(9));
g2_v(17)=(-(getPowerDeriv(y(15),params(7),2)/T(14)));
g2_v(18)=(-((-(T(2)*getPowerDeriv(y(38),params(7),2)))/T(14)));
g2_v(19)=T(17);
g2_v(20)=g2_v(19);
g2_v(21)=(-((-T(10))/T(14)));
g2_v(22)=(-(params(8)*exp(y(23))*1/(1-params(3))));
g2_v(23)=g2_v(22);
g2_v(24)=(-(params(8)*exp(y(23))*(y(18)/(1-params(3))-1)));
g2_v(25)=g2_v(24);
g2_v(26)=(-(params(8)*exp(y(23))*y(13)*1/(1-params(3))));
g2_v(27)=g2_v(26);
g2_v(28)=(-(params(8)*exp(y(23))*exp(y(21))));
g2_v(29)=(-(params(8)*exp(y(23))*exp(y(21))));
g2_v(30)=g2_v(29);
g2_v(31)=(-(T(3)+params(9)*(-(params(8)*exp(y(23))))));
g2_v(32)=1-params(3);
g2_v(33)=g2_v(32);
g2_v(34)=1-params(3);
g2_v(35)=g2_v(34);
g2_v(36)=(-((-1)/(y(17)*y(17))));
g2_v(37)=g2_v(36);
g2_v(38)=(-((-((-y(16))*(y(17)+y(17))))/(y(17)*y(17)*y(17)*y(17))));
g2_v(39)=(-(T(13)*getPowerDeriv(y(18),1-params(10),2)));
g2_v(40)=(-(T(13)*getPowerDeriv(y(18),(-params(10)),2)));
g2_v(41)=(-((T(4)*(-((-y(14))*(y(1)+y(1))))/(y(1)*y(1)*y(1)*y(1))-T(15)*T(15))/(T(4)*T(4))));
g2_v(42)=(-((T(4)*(-1)/(y(1)*y(1))-T(15)*T(16))/(T(4)*T(4))));
g2_v(43)=g2_v(42);
g2_v(44)=(-((-(T(16)*T(16)))/(T(4)*T(4))));
g2_v(45)=(-(100*(-(1/(steady_state(5))*1/(steady_state(5))))/(T(5)*T(5))));
g2_v(46)=(-(100*(-(1/(steady_state(4))*1/(steady_state(4))))/(T(6)*T(6))));
g2_v(47)=(-(100*(-(1/(steady_state(6))*1/(steady_state(6))))/(T(7)*T(7))));
g2_v(48)=(-(100*(-(1/(steady_state(3))*1/(steady_state(3))))/(T(8)*T(8))));
g2 = sparse(g2_i,g2_j,g2_v,23,1849);
end
