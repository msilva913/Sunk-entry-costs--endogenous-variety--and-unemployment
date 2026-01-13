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
    T = BGM.dynamic_g2_tt(T, y, x, params, steady_state, it_);
end
g2_i = zeros(51,1);
g2_j = zeros(51,1);
g2_v = zeros(51,1);

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
g2_i(14)=2;
g2_i(15)=3;
g2_i(16)=3;
g2_i(17)=3;
g2_i(18)=3;
g2_i(19)=3;
g2_i(20)=3;
g2_i(21)=5;
g2_i(22)=5;
g2_i(23)=6;
g2_i(24)=6;
g2_i(25)=6;
g2_i(26)=6;
g2_i(27)=7;
g2_i(28)=7;
g2_i(29)=7;
g2_i(30)=8;
g2_i(31)=8;
g2_i(32)=8;
g2_i(33)=9;
g2_i(34)=9;
g2_i(35)=9;
g2_i(36)=9;
g2_i(37)=9;
g2_i(38)=9;
g2_i(39)=10;
g2_i(40)=10;
g2_i(41)=11;
g2_i(42)=11;
g2_i(43)=11;
g2_i(44)=13;
g2_i(45)=13;
g2_i(46)=14;
g2_i(47)=15;
g2_i(48)=16;
g2_i(49)=17;
g2_i(50)=18;
g2_i(51)=19;
g2_j(1)=82;
g2_j(2)=595;
g2_j(3)=596;
g2_j(4)=621;
g2_j(5)=580;
g2_j(6)=205;
g2_j(7)=597;
g2_j(8)=647;
g2_j(9)=606;
g2_j(10)=206;
g2_j(11)=190;
g2_j(12)=207;
g2_j(13)=632;
g2_j(14)=136;
g2_j(15)=275;
g2_j(16)=375;
g2_j(17)=281;
g2_j(18)=531;
g2_j(19)=385;
g2_j(20)=535;
g2_j(21)=385;
g2_j(22)=535;
g2_j(23)=82;
g2_j(24)=87;
g2_j(25)=212;
g2_j(26)=217;
g2_j(27)=229;
g2_j(28)=529;
g2_j(29)=541;
g2_j(30)=84;
g2_j(31)=134;
g2_j(32)=136;
g2_j(33)=112;
g2_j(34)=187;
g2_j(35)=137;
g2_j(36)=162;
g2_j(37)=218;
g2_j(38)=243;
g2_j(39)=112;
g2_j(40)=187;
g2_j(41)=327;
g2_j(42)=377;
g2_j(43)=379;
g2_j(44)=55;
g2_j(45)=541;
g2_j(46)=325;
g2_j(47)=352;
g2_j(48)=244;
g2_j(49)=271;
g2_j(50)=298;
g2_j(51)=541;
g2_v(1)=T(8);
g2_v(2)=(-((y(25)+y(24))*params(1)*(1-params(3))*getPowerDeriv(y(23),(-params(8)),2)/y(8)));
g2_v(3)=(-(T(7)/y(8)));
g2_v(4)=g2_v(3);
g2_v(5)=(-((-((y(25)+y(24))*T(7)))/(y(8)*y(8))));
g2_v(6)=g2_v(5);
g2_v(7)=(-(T(7)/y(8)));
g2_v(8)=g2_v(7);
g2_v(9)=(-((-T(3))/(y(8)*y(8))));
g2_v(10)=g2_v(9);
g2_v(11)=(-((-((-(T(3)*(y(25)+y(24))))*(y(8)+y(8))))/(y(8)*y(8)*y(8)*y(8))));
g2_v(12)=(-((-T(3))/(y(8)*y(8))));
g2_v(13)=g2_v(12);
g2_v(14)=(-(getPowerDeriv(y(6),1/(params(6)-1),2)));
g2_v(15)=(-y(21));
g2_v(16)=g2_v(15);
g2_v(17)=(-y(15));
g2_v(18)=g2_v(17);
g2_v(19)=(-y(11));
g2_v(20)=g2_v(19);
g2_v(21)=(-(1/T(1)));
g2_v(22)=g2_v(21);
g2_v(23)=(-(T(6)*y(9)*T(8)/params(4)+T(5)*T(5)*T(9)));
g2_v(24)=(-(T(6)*T(4)/params(4)+T(5)*T(2)/params(4)*T(9)));
g2_v(25)=g2_v(24);
g2_v(26)=(-(T(2)/params(4)*T(2)/params(4)*T(9)));
g2_v(27)=(-((-params(7))/(y(21)*y(21))));
g2_v(28)=g2_v(27);
g2_v(29)=(-((-((-(y(9)*params(7)))*(y(21)+y(21))))/(y(21)*y(21)*y(21)*y(21))));
g2_v(30)=(-((-params(6))/(params(6)*y(6)*params(6)*y(6))));
g2_v(31)=g2_v(30);
g2_v(32)=(-((-((-(params(6)*y(4)))*(params(6)*params(6)*y(6)+params(6)*params(6)*y(6))))/(params(6)*y(6)*params(6)*y(6)*params(6)*y(6)*params(6)*y(6))));
g2_v(33)=1;
g2_v(34)=g2_v(33);
g2_v(35)=(-1);
g2_v(36)=g2_v(35);
g2_v(37)=(-1);
g2_v(38)=g2_v(37);
g2_v(39)=(-1);
g2_v(40)=g2_v(39);
g2_v(41)=(-((-1)/(y(15)*y(15))));
g2_v(42)=g2_v(41);
g2_v(43)=(-((-((-y(13))*(y(15)+y(15))))/(y(15)*y(15)*y(15)*y(15))));
g2_v(44)=(-(params(9)*(-1)/(y(3)*y(3))));
g2_v(45)=(-1)/(y(21)*y(21));
g2_v(46)=(-(100*(-1)/(y(13)*y(13))));
g2_v(47)=(-(100*(-1)/(y(14)*y(14))));
g2_v(48)=(-(100*(-1)/(y(10)*y(10))));
g2_v(49)=(-(100*(-1)/(y(11)*y(11))));
g2_v(50)=(-(100*(-1)/(y(12)*y(12))));
g2_v(51)=(-(100*(-1)/(y(21)*y(21))));
g2 = sparse(g2_i,g2_j,g2_v,19,676);
end
