function [residual, g1, g2, g3] = BGM_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [M_.exo_nbr by nperiods] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(19, 1);
mu__ = params(6)/(params(6)-1);
T14 = y(4)^(-params(8));
T21 = params(1)*(1-params(3))*y(23)^(-params(8));
T49 = T14*y(9)/params(4);
T61 = 1/y(9);
T111 = getPowerDeriv(y(4),(-params(8)),1);
T113 = y(9)*T111/params(4);
T114 = getPowerDeriv(T49,params(5),1);
T139 = (-(T21/y(8)));
T147 = T14/params(4);
T153 = (-1)/(y(9)*y(9));
T203 = (-(params(1)*(1-params(3))*getPowerDeriv(y(23),(-params(8)),1)/y(8)));
lhs =T14;
rhs =T21*(y(25)+y(24))/y(8);
residual(1)= lhs-rhs;
lhs =y(15);
rhs =y(6)^(1/(params(6)-1));
residual(2)= lhs-rhs;
lhs =y(4);
rhs =y(15)*y(21)*y(11);
residual(3)= lhs-rhs;
lhs =y(12);
rhs =y(10)-y(11);
residual(4)= lhs-rhs;
lhs =y(9);
rhs =y(15)*y(21)/mu__;
residual(5)= lhs-rhs;
lhs =y(10);
rhs =T49^params(5);
residual(6)= lhs-rhs;
lhs =y(8);
rhs =y(9)*params(7)/y(21);
residual(7)= lhs-rhs;
lhs =y(7);
rhs =y(4)/(params(6)*y(6));
residual(8)= lhs-rhs;
lhs =y(10);
rhs =T61*(y(4)/mu__+y(8)*y(5));
residual(9)= lhs-rhs;
lhs =y(13);
rhs =y(4)+y(8)*y(5);
residual(10)= lhs-rhs;
lhs =y(14);
rhs =y(13)/y(15);
residual(11)= lhs-rhs;
lhs =y(6);
rhs =(1-params(3))*(y(2)+y(1));
residual(12)= lhs-rhs;
lhs =log(y(21));
rhs =params(9)*log(y(3))+x(it_, 1);
residual(13)= lhs-rhs;
lhs =y(16);
rhs =100*log(y(13));
residual(14)= lhs-rhs;
lhs =y(17);
rhs =100*log(y(14));
residual(15)= lhs-rhs;
lhs =y(18);
rhs =100*log(y(10));
residual(16)= lhs-rhs;
lhs =y(19);
rhs =100*log(y(11));
residual(17)= lhs-rhs;
lhs =y(20);
rhs =100*log(y(12));
residual(18)= lhs-rhs;
lhs =y(22);
rhs =log(y(21))*100;
residual(19)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(19, 26);

  %
  % Jacobian matrix
  %

  g1(1,4)=T111;
  g1(1,23)=(-((y(25)+y(24))*params(1)*(1-params(3))*getPowerDeriv(y(23),(-params(8)),1)/y(8)));
  g1(1,24)=T139;
  g1(1,8)=(-((-(T21*(y(25)+y(24))))/(y(8)*y(8))));
  g1(1,25)=T139;
  g1(2,6)=(-(getPowerDeriv(y(6),1/(params(6)-1),1)));
  g1(2,15)=1;
  g1(3,4)=1;
  g1(3,11)=(-(y(15)*y(21)));
  g1(3,15)=(-(y(21)*y(11)));
  g1(3,21)=(-(y(15)*y(11)));
  g1(4,10)=(-1);
  g1(4,11)=1;
  g1(4,12)=1;
  g1(5,9)=1;
  g1(5,15)=(-(y(21)/mu__));
  g1(5,21)=(-(y(15)/mu__));
  g1(6,4)=(-(T113*T114));
  g1(6,9)=(-(T114*T147));
  g1(6,10)=1;
  g1(7,8)=1;
  g1(7,9)=(-(params(7)/y(21)));
  g1(7,21)=(-((-(y(9)*params(7)))/(y(21)*y(21))));
  g1(8,4)=(-(1/(params(6)*y(6))));
  g1(8,6)=(-((-(params(6)*y(4)))/(params(6)*y(6)*params(6)*y(6))));
  g1(8,7)=1;
  g1(9,4)=(-(T61*1/mu__));
  g1(9,5)=(-(y(8)*T61));
  g1(9,8)=(-(T61*y(5)));
  g1(9,9)=(-((y(4)/mu__+y(8)*y(5))*T153));
  g1(9,10)=1;
  g1(10,4)=(-1);
  g1(10,5)=(-y(8));
  g1(10,8)=(-y(5));
  g1(10,13)=1;
  g1(11,13)=(-(1/y(15)));
  g1(11,14)=1;
  g1(11,15)=(-((-y(13))/(y(15)*y(15))));
  g1(12,1)=(-(1-params(3)));
  g1(12,2)=(-(1-params(3)));
  g1(12,6)=1;
  g1(13,3)=(-(params(9)*1/y(3)));
  g1(13,21)=1/y(21);
  g1(13,26)=(-1);
  g1(14,13)=(-(100*1/y(13)));
  g1(14,16)=1;
  g1(15,14)=(-(100*1/y(14)));
  g1(15,17)=1;
  g1(16,10)=(-(100*1/y(10)));
  g1(16,18)=1;
  g1(17,11)=(-(100*1/y(11)));
  g1(17,19)=1;
  g1(18,12)=(-(100*1/y(12)));
  g1(18,20)=1;
  g1(19,21)=(-(100*1/y(21)));
  g1(19,22)=1;
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  v2 = zeros(54,3);
  v2(1,1)=1;
  v2(1,2)=82;
  v2(1,3)=getPowerDeriv(y(4),(-params(8)),2);
  v2(2,1)=1;
  v2(2,2)=595;
  v2(2,3)=(-((y(25)+y(24))*params(1)*(1-params(3))*getPowerDeriv(y(23),(-params(8)),2)/y(8)));
  v2(3,1)=1;
  v2(3,2)=621;
  v2(3,3)=T203;
  v2(4,1)=1;
  v2(4,2)=596;
  v2(4,3)=  v2(3,3);
  v2(5,1)=1;
  v2(5,2)=205;
  v2(5,3)=(-((-((y(25)+y(24))*params(1)*(1-params(3))*getPowerDeriv(y(23),(-params(8)),1)))/(y(8)*y(8))));
  v2(6,1)=1;
  v2(6,2)=580;
  v2(6,3)=  v2(5,3);
  v2(7,1)=1;
  v2(7,2)=206;
  v2(7,3)=(-((-T21)/(y(8)*y(8))));
  v2(8,1)=1;
  v2(8,2)=606;
  v2(8,3)=  v2(7,3);
  v2(9,1)=1;
  v2(9,2)=190;
  v2(9,3)=(-((-((-(T21*(y(25)+y(24))))*(y(8)+y(8))))/(y(8)*y(8)*y(8)*y(8))));
  v2(10,1)=1;
  v2(10,2)=647;
  v2(10,3)=T203;
  v2(11,1)=1;
  v2(11,2)=597;
  v2(11,3)=  v2(10,3);
  v2(12,1)=1;
  v2(12,2)=632;
  v2(12,3)=(-((-T21)/(y(8)*y(8))));
  v2(13,1)=1;
  v2(13,2)=207;
  v2(13,3)=  v2(12,3);
  v2(14,1)=2;
  v2(14,2)=136;
  v2(14,3)=(-(getPowerDeriv(y(6),1/(params(6)-1),2)));
  v2(15,1)=3;
  v2(15,2)=375;
  v2(15,3)=(-y(21));
  v2(16,1)=3;
  v2(16,2)=275;
  v2(16,3)=  v2(15,3);
  v2(17,1)=3;
  v2(17,2)=531;
  v2(17,3)=(-y(15));
  v2(18,1)=3;
  v2(18,2)=281;
  v2(18,3)=  v2(17,3);
  v2(19,1)=3;
  v2(19,2)=535;
  v2(19,3)=(-y(11));
  v2(20,1)=3;
  v2(20,2)=385;
  v2(20,3)=  v2(19,3);
  v2(21,1)=5;
  v2(21,2)=535;
  v2(21,3)=(-(1/mu__));
  v2(22,1)=5;
  v2(22,2)=385;
  v2(22,3)=  v2(21,3);
  v2(23,1)=6;
  v2(23,2)=82;
  v2(23,3)=(-(T114*y(9)*getPowerDeriv(y(4),(-params(8)),2)/params(4)+T113*T113*getPowerDeriv(T49,params(5),2)));
  v2(24,1)=6;
  v2(24,2)=212;
  v2(24,3)=(-(T147*T113*getPowerDeriv(T49,params(5),2)+T114*T111/params(4)));
  v2(25,1)=6;
  v2(25,2)=87;
  v2(25,3)=  v2(24,3);
  v2(26,1)=6;
  v2(26,2)=217;
  v2(26,3)=(-(T147*T147*getPowerDeriv(T49,params(5),2)));
  v2(27,1)=7;
  v2(27,2)=529;
  v2(27,3)=(-((-params(7))/(y(21)*y(21))));
  v2(28,1)=7;
  v2(28,2)=229;
  v2(28,3)=  v2(27,3);
  v2(29,1)=7;
  v2(29,2)=541;
  v2(29,3)=(-((-((-(y(9)*params(7)))*(y(21)+y(21))))/(y(21)*y(21)*y(21)*y(21))));
  v2(30,1)=8;
  v2(30,2)=134;
  v2(30,3)=(-((-params(6))/(params(6)*y(6)*params(6)*y(6))));
  v2(31,1)=8;
  v2(31,2)=84;
  v2(31,3)=  v2(30,3);
  v2(32,1)=8;
  v2(32,2)=136;
  v2(32,3)=(-((-((-(params(6)*y(4)))*(params(6)*params(6)*y(6)+params(6)*params(6)*y(6))))/(params(6)*y(6)*params(6)*y(6)*params(6)*y(6)*params(6)*y(6))));
  v2(33,1)=9;
  v2(33,2)=187;
  v2(33,3)=(-T61);
  v2(34,1)=9;
  v2(34,2)=112;
  v2(34,3)=  v2(33,3);
  v2(35,1)=9;
  v2(35,2)=212;
  v2(35,3)=(-(1/mu__*T153));
  v2(36,1)=9;
  v2(36,2)=87;
  v2(36,3)=  v2(35,3);
  v2(37,1)=9;
  v2(37,2)=213;
  v2(37,3)=(-(y(8)*T153));
  v2(38,1)=9;
  v2(38,2)=113;
  v2(38,3)=  v2(37,3);
  v2(39,1)=9;
  v2(39,2)=216;
  v2(39,3)=(-(y(5)*T153));
  v2(40,1)=9;
  v2(40,2)=191;
  v2(40,3)=  v2(39,3);
  v2(41,1)=9;
  v2(41,2)=217;
  v2(41,3)=(-((y(4)/mu__+y(8)*y(5))*(y(9)+y(9))/(y(9)*y(9)*y(9)*y(9))));
  v2(42,1)=10;
  v2(42,2)=187;
  v2(42,3)=(-1);
  v2(43,1)=10;
  v2(43,2)=112;
  v2(43,3)=  v2(42,3);
  v2(44,1)=11;
  v2(44,2)=377;
  v2(44,3)=(-((-1)/(y(15)*y(15))));
  v2(45,1)=11;
  v2(45,2)=327;
  v2(45,3)=  v2(44,3);
  v2(46,1)=11;
  v2(46,2)=379;
  v2(46,3)=(-((-((-y(13))*(y(15)+y(15))))/(y(15)*y(15)*y(15)*y(15))));
  v2(47,1)=13;
  v2(47,2)=55;
  v2(47,3)=(-(params(9)*(-1)/(y(3)*y(3))));
  v2(48,1)=13;
  v2(48,2)=541;
  v2(48,3)=(-1)/(y(21)*y(21));
  v2(49,1)=14;
  v2(49,2)=325;
  v2(49,3)=(-(100*(-1)/(y(13)*y(13))));
  v2(50,1)=15;
  v2(50,2)=352;
  v2(50,3)=(-(100*(-1)/(y(14)*y(14))));
  v2(51,1)=16;
  v2(51,2)=244;
  v2(51,3)=(-(100*(-1)/(y(10)*y(10))));
  v2(52,1)=17;
  v2(52,2)=271;
  v2(52,3)=(-(100*(-1)/(y(11)*y(11))));
  v2(53,1)=18;
  v2(53,2)=298;
  v2(53,3)=(-(100*(-1)/(y(12)*y(12))));
  v2(54,1)=19;
  v2(54,2)=541;
  v2(54,3)=(-(100*(-1)/(y(21)*y(21))));
  g2 = sparse(v2(:,1),v2(:,2),v2(:,3),19,676);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],19,17576);
end
end
