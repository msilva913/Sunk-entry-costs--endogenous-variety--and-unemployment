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
    T = BGM_remove_variety_effects.dynamic_g2_tt(T, y, x, params, steady_state, it_);
end
g2_i = zeros(41,1);
g2_j = zeros(41,1);
g2_v = zeros(41,1);

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
g2_i(14)=3;
g2_i(15)=3;
g2_i(16)=6;
g2_i(17)=6;
g2_i(18)=6;
g2_i(19)=6;
g2_i(20)=7;
g2_i(21)=7;
g2_i(22)=7;
g2_i(23)=8;
g2_i(24)=8;
g2_i(25)=8;
g2_i(26)=9;
g2_i(27)=9;
g2_i(28)=9;
g2_i(29)=9;
g2_i(30)=9;
g2_i(31)=9;
g2_i(32)=10;
g2_i(33)=10;
g2_i(34)=13;
g2_i(35)=13;
g2_i(36)=14;
g2_i(37)=15;
g2_i(38)=16;
g2_i(39)=17;
g2_i(40)=18;
g2_i(41)=19;
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
g2_j(14)=281;
g2_j(15)=531;
g2_j(16)=82;
g2_j(17)=87;
g2_j(18)=212;
g2_j(19)=217;
g2_j(20)=229;
g2_j(21)=529;
g2_j(22)=541;
g2_j(23)=84;
g2_j(24)=134;
g2_j(25)=136;
g2_j(26)=112;
g2_j(27)=187;
g2_j(28)=137;
g2_j(29)=162;
g2_j(30)=218;
g2_j(31)=243;
g2_j(32)=112;
g2_j(33)=187;
g2_j(34)=55;
g2_j(35)=541;
g2_j(36)=325;
g2_j(37)=352;
g2_j(38)=244;
g2_j(39)=271;
g2_j(40)=298;
g2_j(41)=541;
g2_v(1)=T(7);
g2_v(2)=(-((y(25)+y(24))*params(1)*(1-params(3))*getPowerDeriv(y(23),(-params(8)),2)/y(8)));
g2_v(3)=(-(T(6)/y(8)));
g2_v(4)=g2_v(3);
g2_v(5)=(-((-((y(25)+y(24))*T(6)))/(y(8)*y(8))));
g2_v(6)=g2_v(5);
g2_v(7)=(-(T(6)/y(8)));
g2_v(8)=g2_v(7);
g2_v(9)=(-((-T(2))/(y(8)*y(8))));
g2_v(10)=g2_v(9);
g2_v(11)=(-((-((-(T(2)*(y(25)+y(24))))*(y(8)+y(8))))/(y(8)*y(8)*y(8)*y(8))));
g2_v(12)=(-((-T(2))/(y(8)*y(8))));
g2_v(13)=g2_v(12);
g2_v(14)=(-1.0);
g2_v(15)=g2_v(14);
g2_v(16)=(-(T(5)*y(9)*T(7)/params(4)+T(4)*T(4)*T(8)));
g2_v(17)=(-(T(5)*T(3)/params(4)+T(4)*T(1)/params(4)*T(8)));
g2_v(18)=g2_v(17);
g2_v(19)=(-(T(1)/params(4)*T(1)/params(4)*T(8)));
g2_v(20)=(-((-params(7))/(y(21)*y(21))));
g2_v(21)=g2_v(20);
g2_v(22)=(-((-((-(y(9)*params(7)))*(y(21)+y(21))))/(y(21)*y(21)*y(21)*y(21))));
g2_v(23)=(-((-params(6))/(params(6)*y(6)*params(6)*y(6))));
g2_v(24)=g2_v(23);
g2_v(25)=(-((-((-(params(6)*y(4)))*(params(6)*params(6)*y(6)+params(6)*params(6)*y(6))))/(params(6)*y(6)*params(6)*y(6)*params(6)*y(6)*params(6)*y(6))));
g2_v(26)=1;
g2_v(27)=g2_v(26);
g2_v(28)=(-1);
g2_v(29)=g2_v(28);
g2_v(30)=(-1);
g2_v(31)=g2_v(30);
g2_v(32)=(-1);
g2_v(33)=g2_v(32);
g2_v(34)=(-(params(9)*(-1)/(y(3)*y(3))));
g2_v(35)=(-1)/(y(21)*y(21));
g2_v(36)=(-(100*(-1)/(y(13)*y(13))));
g2_v(37)=(-(100*(-1)/(y(14)*y(14))));
g2_v(38)=(-(100*(-1)/(y(10)*y(10))));
g2_v(39)=(-(100*(-1)/(y(11)*y(11))));
g2_v(40)=(-(100*(-1)/(y(12)*y(12))));
g2_v(41)=(-(100*(-1)/(y(21)*y(21))));
g2 = sparse(g2_i,g2_j,g2_v,19,676);
end
