function [residual, g1, g2, g3] = model_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(65, 1);
mu__ = params(13)/(params(13)-1);
w_ss__ = (steady_state(13));
T48 = params(16)*(1-y(69))*y(106)/y(44);
T80 = y(111)-y(112)-y(109)+(1-y(113))*(params(17)+y(109)/y(105));
T332 = 0.3333333333333333*(y(41)*(1-y(69))+(1-y(11))*y(1)+(1-y(22))*y(23));
T459 = (-(params(16)*(1-y(69))*y(106)))/(y(44)*y(44));
T477 = (-(0.3333333333333333/(0.3333333333333333*(y(47)+y(3)+y(24)))));
T517 = (-(0.3333333333333333/(0.3333333333333333*(y(5)+y(54)+y(19)))));
T521 = (-(0.3333333333333333/(0.3333333333333333*(y(55)+y(6)+y(20)))));
T527 = (-(0.3333333333333333/(0.3333333333333333*(y(56)+y(7)+y(21)))));
T547 = (-(0.3333333333333333/(0.3333333333333333*(y(66)+y(8)+y(25)))));
T574 = (-(1/(y(69)+y(11)+y(22))));
lhs =y(40);
rhs =y(55)/y(5);
residual(1)= lhs-rhs;
lhs =y(41);
rhs =params(19)*y(40)^(1-params(11));
residual(2)= lhs-rhs;
lhs =y(42);
rhs =params(19)*y(40)^(-params(11));
residual(3)= lhs-rhs;
lhs =y(43);
rhs =y(4)^params(14);
residual(4)= lhs-rhs;
lhs =y(44);
rhs =y(61)^(-params(1));
residual(5)= lhs-rhs;
lhs =y(45);
rhs =T48*(y(107)+y(108));
residual(6)= lhs-rhs;
lhs =y(45)*y(47);
rhs =y(51)*y(58);
residual(7)= lhs-rhs;
lhs =y(47)*params(18);
rhs =y(58)*y(67);
residual(8)= lhs-rhs;
lhs =params(17)+y(48)/y(42);
rhs =T48*T80;
residual(9)= lhs-rhs;
lhs =y(48);
rhs =y(49)-T48*y(110);
residual(10)= lhs-rhs;
lhs =y(49);
rhs =params(20)*(y(50)/params(12))^params(3);
residual(11)= lhs-rhs;
lhs =y(51);
rhs =y(43)*y(67)/mu__;
residual(12)= lhs-rhs;
lhs =y(52);
rhs =params(21)*(y(51)-y(48)+y(41)*(params(17)+y(48)/y(42)))+(1-params(21))*params(2)*w_ss__;
residual(13)= lhs-rhs;
lhs =y(53);
rhs =(1-y(69))*(y(4)+y(47));
residual(14)= lhs-rhs;
lhs =y(54);
rhs =y(5)*(1-y(41)*(1-y(69)))+y(56)*(1-y(5));
residual(15)= lhs-rhs;
lhs =y(55);
rhs =y(50)+(1-y(11))*((1-y(2))*y(6)+(1-y(5))*y(10));
residual(16)= lhs-rhs;
lhs =y(56);
rhs =1-(1-y(69))*(1-y(68));
residual(17)= lhs-rhs;
lhs =y(57);
rhs =1-y(5);
residual(18)= lhs-rhs;
lhs =y(57);
rhs =y(58)+y(59);
residual(19)= lhs-rhs;
lhs =y(62);
rhs =y(43)*y(67)*y(59);
residual(20)= lhs-rhs;
lhs =y(60);
rhs =y(49)*y(50)/(1+params(3))+y(42)*y(55)*params(17);
residual(21)= lhs-rhs;
lhs =y(62);
rhs =y(61)+y(60);
residual(22)= lhs-rhs;
lhs =y(63);
rhs =y(45)*y(47)+y(62);
residual(23)= lhs-rhs;
lhs =y(63);
rhs =y(51)*y(57)+y(4)*y(46);
residual(24)= lhs-rhs;
lhs =y(64);
rhs =y(63)-y(60);
residual(25)= lhs-rhs;
lhs =y(65);
rhs =y(52)*y(57)/y(64);
residual(26)= lhs-rhs;
lhs =y(66);
rhs =y(64)/(y(43)*y(57));
residual(27)= lhs-rhs;
lhs =log(y(67));
rhs =params(5)*log(y(9))+(1-params(5))*log((steady_state(28)))-params(8)*x(it_, 1);
residual(28)= lhs-rhs;
lhs =log(y(68));
rhs =params(6)*log(y(10))+(1-params(6))*log((steady_state(29)))+params(9)*x(it_, 2);
residual(29)= lhs-rhs;
lhs =log(y(69));
rhs =params(7)*log(y(11))+(1-params(7))*log((steady_state(30)))+params(10)*x(it_, 3);
residual(30)= lhs-rhs;
lhs =y(70);
rhs =log(0.3333333333333333*(y(5)+y(54)+y(19)));
residual(31)= lhs-rhs;
lhs =y(71);
rhs =log(0.3333333333333333*(y(55)+y(6)+y(20)));
residual(32)= lhs-rhs;
lhs =y(72);
rhs =log(0.3333333333333333*(y(56)+y(7)+y(21)));
residual(33)= lhs-rhs;
lhs =y(73);
rhs =log(T332);
residual(34)= lhs-rhs;
lhs =y(74);
rhs =log(y(69)+y(11)+y(22));
residual(35)= lhs-rhs;
lhs =y(75);
rhs =log(0.3333333333333333*(y(47)+y(3)+y(24)));
residual(36)= lhs-rhs;
lhs =y(76);
rhs =log(0.3333333333333333*(y(66)+y(8)+y(25)));
residual(37)= lhs-rhs;
lhs =y(77);
rhs =y(27);
residual(38)= lhs-rhs;
lhs =y(78);
rhs =y(29);
residual(39)= lhs-rhs;
lhs =y(79);
rhs =y(31);
residual(40)= lhs-rhs;
lhs =y(80);
rhs =y(33);
residual(41)= lhs-rhs;
lhs =y(81);
rhs =y(35);
residual(42)= lhs-rhs;
lhs =y(82);
rhs =y(37);
residual(43)= lhs-rhs;
lhs =y(83);
rhs =y(39);
residual(44)= lhs-rhs;
lhs =y(84);
rhs =y(5);
residual(45)= lhs-rhs;
lhs =y(85);
rhs =y(6);
residual(46)= lhs-rhs;
lhs =y(86);
rhs =y(7);
residual(47)= lhs-rhs;
lhs =y(87);
rhs =y(11);
residual(48)= lhs-rhs;
lhs =y(88);
rhs =y(1);
residual(49)= lhs-rhs;
lhs =y(89);
rhs =y(3);
residual(50)= lhs-rhs;
lhs =y(90);
rhs =y(8);
residual(51)= lhs-rhs;
lhs =y(91);
rhs =y(12);
residual(52)= lhs-rhs;
lhs =y(92);
rhs =y(26);
residual(53)= lhs-rhs;
lhs =y(93);
rhs =y(13);
residual(54)= lhs-rhs;
lhs =y(94);
rhs =y(28);
residual(55)= lhs-rhs;
lhs =y(95);
rhs =y(14);
residual(56)= lhs-rhs;
lhs =y(96);
rhs =y(30);
residual(57)= lhs-rhs;
lhs =y(97);
rhs =y(15);
residual(58)= lhs-rhs;
lhs =y(98);
rhs =y(32);
residual(59)= lhs-rhs;
lhs =y(99);
rhs =y(16);
residual(60)= lhs-rhs;
lhs =y(100);
rhs =y(34);
residual(61)= lhs-rhs;
lhs =y(101);
rhs =y(17);
residual(62)= lhs-rhs;
lhs =y(102);
rhs =y(36);
residual(63)= lhs-rhs;
lhs =y(103);
rhs =y(18);
residual(64)= lhs-rhs;
lhs =y(104);
rhs =y(38);
residual(65)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(65, 116);

  %
  % Jacobian matrix
  %

  g1(1,40)=1;
  g1(1,5)=(-((-y(55))/(y(5)*y(5))));
  g1(1,55)=(-(1/y(5)));
  g1(2,40)=(-(params(19)*getPowerDeriv(y(40),1-params(11),1)));
  g1(2,41)=1;
  g1(3,40)=(-(params(19)*getPowerDeriv(y(40),(-params(11)),1)));
  g1(3,42)=1;
  g1(4,43)=1;
  g1(4,4)=(-(getPowerDeriv(y(4),params(14),1)));
  g1(5,44)=1;
  g1(5,61)=(-(getPowerDeriv(y(61),(-params(1)),1)));
  g1(6,44)=(-((y(107)+y(108))*T459));
  g1(6,106)=(-((y(107)+y(108))*params(16)*(1-y(69))/y(44)));
  g1(6,45)=1;
  g1(6,107)=(-T48);
  g1(6,108)=(-T48);
  g1(6,69)=(-((y(107)+y(108))*y(106)*(-params(16))/y(44)));
  g1(7,45)=y(47);
  g1(7,47)=y(45);
  g1(7,51)=(-y(58));
  g1(7,58)=(-y(51));
  g1(8,47)=params(18);
  g1(8,58)=(-y(67));
  g1(8,67)=(-y(58));
  g1(9,42)=(-y(48))/(y(42)*y(42));
  g1(9,105)=(-(T48*(1-y(113))*(-y(109))/(y(105)*y(105))));
  g1(9,44)=(-(T80*T459));
  g1(9,106)=(-(T80*params(16)*(1-y(69))/y(44)));
  g1(9,48)=1/y(42);
  g1(9,109)=(-(T48*((-1)+(1-y(113))*1/y(105))));
  g1(9,111)=(-T48);
  g1(9,112)=T48;
  g1(9,113)=(-(T48*(-(params(17)+y(109)/y(105)))));
  g1(9,69)=(-(T80*y(106)*(-params(16))/y(44)));
  g1(10,44)=y(110)*T459;
  g1(10,106)=y(110)*params(16)*(1-y(69))/y(44);
  g1(10,48)=1;
  g1(10,49)=(-1);
  g1(10,110)=T48;
  g1(10,69)=y(110)*y(106)*(-params(16))/y(44);
  g1(11,49)=1;
  g1(11,50)=(-(params(20)*1/params(12)*getPowerDeriv(y(50)/params(12),params(3),1)));
  g1(12,43)=(-(y(67)/mu__));
  g1(12,51)=1;
  g1(12,67)=(-(y(43)/mu__));
  g1(13,41)=(-((params(17)+y(48)/y(42))*params(21)));
  g1(13,42)=(-(params(21)*y(41)*(-y(48))/(y(42)*y(42))));
  g1(13,48)=(-(params(21)*((-1)+y(41)*1/y(42))));
  g1(13,51)=(-params(21));
  g1(13,52)=1;
  g1(14,47)=(-(1-y(69)));
  g1(14,4)=(-(1-y(69)));
  g1(14,53)=1;
  g1(14,69)=y(4)+y(47);
  g1(15,41)=(-(y(5)*(-(1-y(69)))));
  g1(15,5)=(-(1-y(41)*(1-y(69))-y(56)));
  g1(15,54)=1;
  g1(15,56)=(-(1-y(5)));
  g1(15,69)=(-(y(5)*y(41)));
  g1(16,2)=(-((1-y(11))*(-y(6))));
  g1(16,50)=(-1);
  g1(16,5)=(-((1-y(11))*(-y(10))));
  g1(16,6)=(-((1-y(11))*(1-y(2))));
  g1(16,55)=1;
  g1(16,10)=(-((1-y(5))*(1-y(11))));
  g1(16,11)=(1-y(2))*y(6)+(1-y(5))*y(10);
  g1(17,56)=1;
  g1(17,68)=(-(1-y(69)));
  g1(17,69)=(-(1-y(68)));
  g1(18,5)=1;
  g1(18,57)=1;
  g1(19,57)=1;
  g1(19,58)=(-1);
  g1(19,59)=(-1);
  g1(20,43)=(-(y(67)*y(59)));
  g1(20,59)=(-(y(43)*y(67)));
  g1(20,62)=1;
  g1(20,67)=(-(y(43)*y(59)));
  g1(21,42)=(-(y(55)*params(17)));
  g1(21,49)=(-(y(50)/(1+params(3))));
  g1(21,50)=(-(y(49)/(1+params(3))));
  g1(21,55)=(-(y(42)*params(17)));
  g1(21,60)=1;
  g1(22,60)=(-1);
  g1(22,61)=(-1);
  g1(22,62)=1;
  g1(23,45)=(-y(47));
  g1(23,47)=(-y(45));
  g1(23,62)=(-1);
  g1(23,63)=1;
  g1(24,46)=(-y(4));
  g1(24,51)=(-y(57));
  g1(24,4)=(-y(46));
  g1(24,57)=(-y(51));
  g1(24,63)=1;
  g1(25,60)=1;
  g1(25,63)=(-1);
  g1(25,64)=1;
  g1(26,52)=(-(y(57)/y(64)));
  g1(26,57)=(-(y(52)/y(64)));
  g1(26,64)=(-((-(y(52)*y(57)))/(y(64)*y(64))));
  g1(26,65)=1;
  g1(27,43)=(-((-(y(57)*y(64)))/(y(43)*y(57)*y(43)*y(57))));
  g1(27,57)=(-((-(y(43)*y(64)))/(y(43)*y(57)*y(43)*y(57))));
  g1(27,64)=(-(1/(y(43)*y(57))));
  g1(27,66)=1;
  g1(28,9)=(-(params(5)*1/y(9)));
  g1(28,67)=1/y(67);
  g1(28,114)=params(8);
  g1(29,10)=(-(params(6)*1/y(10)));
  g1(29,68)=1/y(68);
  g1(29,115)=(-params(9));
  g1(30,11)=(-(params(7)*1/y(11)));
  g1(30,69)=1/y(69);
  g1(30,116)=(-params(10));
  g1(31,5)=T517;
  g1(31,54)=T517;
  g1(31,70)=1;
  g1(31,19)=T517;
  g1(32,6)=T521;
  g1(32,55)=T521;
  g1(32,71)=1;
  g1(32,20)=T521;
  g1(33,7)=T527;
  g1(33,56)=T527;
  g1(33,72)=1;
  g1(33,21)=T527;
  g1(34,1)=(-((1-y(11))*0.3333333333333333/T332));
  g1(34,41)=(-((1-y(69))*0.3333333333333333/T332));
  g1(34,11)=(-(0.3333333333333333*(-y(1))/T332));
  g1(34,69)=(-(0.3333333333333333*(-y(41))/T332));
  g1(34,73)=1;
  g1(34,22)=(-(0.3333333333333333*(-y(23))/T332));
  g1(34,23)=(-(0.3333333333333333*(1-y(22))/T332));
  g1(35,11)=T574;
  g1(35,69)=T574;
  g1(35,74)=1;
  g1(35,22)=T574;
  g1(36,3)=T477;
  g1(36,47)=T477;
  g1(36,75)=1;
  g1(36,24)=T477;
  g1(37,8)=T547;
  g1(37,66)=T547;
  g1(37,76)=1;
  g1(37,25)=T547;
  g1(38,77)=1;
  g1(38,27)=(-1);
  g1(39,78)=1;
  g1(39,29)=(-1);
  g1(40,79)=1;
  g1(40,31)=(-1);
  g1(41,80)=1;
  g1(41,33)=(-1);
  g1(42,81)=1;
  g1(42,35)=(-1);
  g1(43,82)=1;
  g1(43,37)=(-1);
  g1(44,83)=1;
  g1(44,39)=(-1);
  g1(45,5)=(-1);
  g1(45,84)=1;
  g1(46,6)=(-1);
  g1(46,85)=1;
  g1(47,7)=(-1);
  g1(47,86)=1;
  g1(48,11)=(-1);
  g1(48,87)=1;
  g1(49,1)=(-1);
  g1(49,88)=1;
  g1(50,3)=(-1);
  g1(50,89)=1;
  g1(51,8)=(-1);
  g1(51,90)=1;
  g1(52,12)=(-1);
  g1(52,91)=1;
  g1(53,26)=(-1);
  g1(53,92)=1;
  g1(54,13)=(-1);
  g1(54,93)=1;
  g1(55,28)=(-1);
  g1(55,94)=1;
  g1(56,14)=(-1);
  g1(56,95)=1;
  g1(57,30)=(-1);
  g1(57,96)=1;
  g1(58,15)=(-1);
  g1(58,97)=1;
  g1(59,32)=(-1);
  g1(59,98)=1;
  g1(60,16)=(-1);
  g1(60,99)=1;
  g1(61,34)=(-1);
  g1(61,100)=1;
  g1(62,17)=(-1);
  g1(62,101)=1;
  g1(63,36)=(-1);
  g1(63,102)=1;
  g1(64,18)=(-1);
  g1(64,103)=1;
  g1(65,38)=(-1);
  g1(65,104)=1;
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],65,13456);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],65,1560896);
end
end
