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

residual = zeros(63, 1);
T47 = 1-exp(y(67));
T70 = T47*params(12)*y(103)/y(43);
T83 = y(108)-y(109)-y(106)+(1-exp(y(66)))*(params(18)+y(106)/y(102));
T93 = y(49)/params(19);
T105 = y(39)/T47;
T459 = (-(1/(y(46)+y(3)+y(24))));
T493 = (-(0.3333333333333333/(0.3333333333333333*(y(53)+y(5)+y(19)))));
T501 = (-(0.3333333333333333/(0.3333333333333333*(y(54)+y(6)+y(20)))));
T524 = (-(0.3333333333333333/(0.3333333333333333*(y(63)+y(7)+y(21)))));
T526 = (-(1/(y(64)+y(8)+y(22))));
lhs =y(39);
rhs =y(54)/y(53);
residual(1)= lhs-rhs;
lhs =y(40);
rhs =params(16)*y(39)^(1-params(13));
residual(2)= lhs-rhs;
lhs =y(41);
rhs =params(16)*y(39)^(-params(13));
residual(3)= lhs-rhs;
lhs =y(42);
rhs =y(52)^(1/(params(5)-1));
residual(4)= lhs-rhs;
lhs =y(43);
rhs =y(59)^(-params(1));
residual(5)= lhs-rhs;
lhs =y(44);
rhs =(params(5)-1)*y(42)*params(17)/params(5);
residual(6)= lhs-rhs;
lhs =y(44);
rhs =params(12)*T47*y(103)/y(43)*(y(104)+y(105));
residual(7)= lhs-rhs;
lhs =y(46);
rhs =y(56)*exp(y(65))/params(17);
residual(8)= lhs-rhs;
lhs =params(18)+y(47)/y(41);
rhs =T70*T83;
residual(9)= lhs-rhs;
lhs =y(47);
rhs =y(48)-T70*y(107);
residual(10)= lhs-rhs;
lhs =y(48);
rhs =T93^params(3);
residual(11)= lhs-rhs;
lhs =y(50);
rhs =(params(5)-1)*y(42)*exp(y(65))/params(5);
residual(12)= lhs-rhs;
lhs =y(51);
rhs =params(20)*(y(50)-y(47)+T105*(y(47)+y(41)*params(18)))+(1-params(20))*params(2)*(steady_state(13));
residual(13)= lhs-rhs;
lhs =y(52);
rhs =T47*(y(4)+y(3));
residual(14)= lhs-rhs;
lhs =y(53);
rhs =(1-T47*y(1))*y(5)+y(64)*(1-y(5));
residual(15)= lhs-rhs;
lhs =y(54);
rhs =y(49)+T47*((1-y(2))*y(6)+exp(y(66))*(1-y(5)));
residual(16)= lhs-rhs;
lhs =y(55);
rhs =1-y(53);
residual(17)= lhs-rhs;
lhs =y(55);
rhs =y(56)+y(57);
residual(18)= lhs-rhs;
lhs =y(60);
rhs =y(42)*exp(y(65))*y(57);
residual(19)= lhs-rhs;
lhs =y(58);
rhs =params(19)/(1+params(3))*T93^(1+params(3));
residual(20)= lhs-rhs;
lhs =y(60);
rhs =y(59)+y(58)+y(41)*y(54)*params(18);
residual(21)= lhs-rhs;
lhs =y(61);
rhs =y(60)+y(44)*y(46);
residual(22)= lhs-rhs;
lhs =y(61);
rhs =y(50)*y(55)+y(52)*y(45);
residual(23)= lhs-rhs;
lhs =y(62);
rhs =y(51)*y(55)/y(61);
residual(24)= lhs-rhs;
lhs =y(63);
rhs =y(61)/(y(42)*y(55));
residual(25)= lhs-rhs;
lhs =y(64);
rhs =1-T47*(1-exp(y(66)));
residual(26)= lhs-rhs;
lhs =y(65);
rhs =params(6)*y(9)+(1-params(6))*(steady_state(27))-params(9)*x(it_, 1);
residual(27)= lhs-rhs;
lhs =y(66);
rhs =params(7)*y(10)+(1-params(7))*(steady_state(28))+params(10)*x(it_, 2);
residual(28)= lhs-rhs;
lhs =y(67);
rhs =params(8)*y(11)+(1-params(8))*(steady_state(29))+params(11)*x(it_, 3);
residual(29)= lhs-rhs;
lhs =y(68);
rhs =log(0.3333333333333333*(y(53)+y(5)+y(19)));
residual(30)= lhs-rhs;
lhs =y(69);
rhs =log(0.3333333333333333*(y(54)+y(6)+y(20)));
residual(31)= lhs-rhs;
lhs =y(70);
rhs =y(69)-y(68);
residual(32)= lhs-rhs;
lhs =y(71);
rhs =log(0.3333333333333333*(y(63)+y(7)+y(21)));
residual(33)= lhs-rhs;
lhs =y(72);
rhs =log(y(64)+y(8)+y(22));
residual(34)= lhs-rhs;
lhs =y(73);
rhs =log(exp(y(67))+exp(y(11))+exp(y(23)));
residual(35)= lhs-rhs;
lhs =y(74);
rhs =log(y(46)+y(3)+y(24));
residual(36)= lhs-rhs;
lhs =y(75);
rhs =y(26);
residual(37)= lhs-rhs;
lhs =y(76);
rhs =y(28);
residual(38)= lhs-rhs;
lhs =y(77);
rhs =y(30);
residual(39)= lhs-rhs;
lhs =y(78);
rhs =y(32);
residual(40)= lhs-rhs;
lhs =y(79);
rhs =y(34);
residual(41)= lhs-rhs;
lhs =y(80);
rhs =y(36);
residual(42)= lhs-rhs;
lhs =y(81);
rhs =y(38);
residual(43)= lhs-rhs;
lhs =y(82);
rhs =y(5);
residual(44)= lhs-rhs;
lhs =y(83);
rhs =y(6);
residual(45)= lhs-rhs;
lhs =y(84);
rhs =y(7);
residual(46)= lhs-rhs;
lhs =y(85);
rhs =y(8);
residual(47)= lhs-rhs;
lhs =y(86);
rhs =y(11);
residual(48)= lhs-rhs;
lhs =y(87);
rhs =y(3);
residual(49)= lhs-rhs;
lhs =y(88);
rhs =y(12);
residual(50)= lhs-rhs;
lhs =y(89);
rhs =y(25);
residual(51)= lhs-rhs;
lhs =y(90);
rhs =y(13);
residual(52)= lhs-rhs;
lhs =y(91);
rhs =y(27);
residual(53)= lhs-rhs;
lhs =y(92);
rhs =y(14);
residual(54)= lhs-rhs;
lhs =y(93);
rhs =y(29);
residual(55)= lhs-rhs;
lhs =y(94);
rhs =y(15);
residual(56)= lhs-rhs;
lhs =y(95);
rhs =y(31);
residual(57)= lhs-rhs;
lhs =y(96);
rhs =y(16);
residual(58)= lhs-rhs;
lhs =y(97);
rhs =y(33);
residual(59)= lhs-rhs;
lhs =y(98);
rhs =y(17);
residual(60)= lhs-rhs;
lhs =y(99);
rhs =y(35);
residual(61)= lhs-rhs;
lhs =y(100);
rhs =y(18);
residual(62)= lhs-rhs;
lhs =y(101);
rhs =y(37);
residual(63)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(63, 112);

  %
  % Jacobian matrix
  %

  g1(1,39)=1;
  g1(1,53)=(-((-y(54))/(y(53)*y(53))));
  g1(1,54)=(-(1/y(53)));
  g1(2,39)=(-(params(16)*getPowerDeriv(y(39),1-params(13),1)));
  g1(2,40)=1;
  g1(3,39)=(-(params(16)*getPowerDeriv(y(39),(-params(13)),1)));
  g1(3,41)=1;
  g1(4,42)=1;
  g1(4,52)=(-(getPowerDeriv(y(52),1/(params(5)-1),1)));
  g1(5,43)=1;
  g1(5,59)=(-(getPowerDeriv(y(59),(-params(1)),1)));
  g1(6,42)=(-((params(5)-1)*params(17)/params(5)));
  g1(6,44)=1;
  g1(7,43)=(-((y(104)+y(105))*(-(params(12)*T47*y(103)))/(y(43)*y(43))));
  g1(7,103)=(-((y(104)+y(105))*params(12)*T47/y(43)));
  g1(7,44)=1;
  g1(7,104)=(-(params(12)*T47*y(103)/y(43)));
  g1(7,105)=(-(params(12)*T47*y(103)/y(43)));
  g1(7,67)=(-((y(104)+y(105))*y(103)*params(12)*(-exp(y(67)))/y(43)));
  g1(8,46)=1;
  g1(8,56)=(-(exp(y(65))/params(17)));
  g1(8,65)=(-(y(56)*exp(y(65))/params(17)));
  g1(9,41)=(-y(47))/(y(41)*y(41));
  g1(9,102)=(-(T70*(1-exp(y(66)))*(-y(106))/(y(102)*y(102))));
  g1(9,43)=(-(T83*T47*(-(params(12)*y(103)))/(y(43)*y(43))));
  g1(9,103)=(-(T83*T47*params(12)/y(43)));
  g1(9,47)=1/y(41);
  g1(9,106)=(-(T70*((-1)+(1-exp(y(66)))*1/y(102))));
  g1(9,108)=(-T70);
  g1(9,109)=T70;
  g1(9,66)=(-(T70*(params(18)+y(106)/y(102))*(-exp(y(66)))));
  g1(9,67)=(-(T83*params(12)*y(103)/y(43)*(-exp(y(67)))));
  g1(10,43)=y(107)*T47*(-(params(12)*y(103)))/(y(43)*y(43));
  g1(10,103)=y(107)*T47*params(12)/y(43);
  g1(10,47)=1;
  g1(10,48)=(-1);
  g1(10,107)=T70;
  g1(10,67)=y(107)*params(12)*y(103)/y(43)*(-exp(y(67)));
  g1(11,48)=1;
  g1(11,49)=(-(1/params(19)*getPowerDeriv(T93,params(3),1)));
  g1(12,42)=(-((params(5)-1)*exp(y(65))/params(5)));
  g1(12,50)=1;
  g1(12,65)=(-((params(5)-1)*y(42)*exp(y(65))/params(5)));
  g1(13,39)=(-(params(20)*(y(47)+y(41)*params(18))*1/T47));
  g1(13,41)=(-(params(20)*params(18)*T105));
  g1(13,47)=(-(params(20)*((-1)+T105)));
  g1(13,50)=(-params(20));
  g1(13,51)=1;
  g1(13,67)=(-(params(20)*(y(47)+y(41)*params(18))*(-(y(39)*(-exp(y(67)))))/(T47*T47)));
  g1(14,3)=(-T47);
  g1(14,4)=(-T47);
  g1(14,52)=1;
  g1(14,67)=(-((y(4)+y(3))*(-exp(y(67)))));
  g1(15,1)=(-(y(5)*(-T47)));
  g1(15,5)=(-(1-T47*y(1)-y(64)));
  g1(15,53)=1;
  g1(15,64)=(-(1-y(5)));
  g1(15,67)=(-(y(5)*(-(y(1)*(-exp(y(67)))))));
  g1(16,2)=(-(T47*(-y(6))));
  g1(16,49)=(-1);
  g1(16,5)=(-(T47*(-exp(y(66)))));
  g1(16,6)=(-(T47*(1-y(2))));
  g1(16,54)=1;
  g1(16,66)=(-(T47*exp(y(66))*(1-y(5))));
  g1(16,67)=(-(((1-y(2))*y(6)+exp(y(66))*(1-y(5)))*(-exp(y(67)))));
  g1(17,53)=1;
  g1(17,55)=1;
  g1(18,55)=1;
  g1(18,56)=(-1);
  g1(18,57)=(-1);
  g1(19,42)=(-(exp(y(65))*y(57)));
  g1(19,57)=(-(y(42)*exp(y(65))));
  g1(19,60)=1;
  g1(19,65)=(-(y(42)*exp(y(65))*y(57)));
  g1(20,49)=(-(params(19)/(1+params(3))*1/params(19)*getPowerDeriv(T93,1+params(3),1)));
  g1(20,58)=1;
  g1(21,41)=(-(y(54)*params(18)));
  g1(21,54)=(-(y(41)*params(18)));
  g1(21,58)=(-1);
  g1(21,59)=(-1);
  g1(21,60)=1;
  g1(22,44)=(-y(46));
  g1(22,46)=(-y(44));
  g1(22,60)=(-1);
  g1(22,61)=1;
  g1(23,45)=(-y(52));
  g1(23,50)=(-y(55));
  g1(23,52)=(-y(45));
  g1(23,55)=(-y(50));
  g1(23,61)=1;
  g1(24,51)=(-(y(55)/y(61)));
  g1(24,55)=(-(y(51)/y(61)));
  g1(24,61)=(-((-(y(51)*y(55)))/(y(61)*y(61))));
  g1(24,62)=1;
  g1(25,42)=(-((-(y(55)*y(61)))/(y(42)*y(55)*y(42)*y(55))));
  g1(25,55)=(-((-(y(42)*y(61)))/(y(42)*y(55)*y(42)*y(55))));
  g1(25,61)=(-(1/(y(42)*y(55))));
  g1(25,63)=1;
  g1(26,64)=1;
  g1(26,66)=T47*(-exp(y(66)));
  g1(26,67)=(1-exp(y(66)))*(-exp(y(67)));
  g1(27,9)=(-params(6));
  g1(27,65)=1;
  g1(27,110)=params(9);
  g1(28,10)=(-params(7));
  g1(28,66)=1;
  g1(28,111)=(-params(10));
  g1(29,11)=(-params(8));
  g1(29,67)=1;
  g1(29,112)=(-params(11));
  g1(30,5)=T493;
  g1(30,53)=T493;
  g1(30,68)=1;
  g1(30,19)=T493;
  g1(31,6)=T501;
  g1(31,54)=T501;
  g1(31,69)=1;
  g1(31,20)=T501;
  g1(32,68)=1;
  g1(32,69)=(-1);
  g1(32,70)=1;
  g1(33,7)=T524;
  g1(33,63)=T524;
  g1(33,71)=1;
  g1(33,21)=T524;
  g1(34,8)=T526;
  g1(34,64)=T526;
  g1(34,72)=1;
  g1(34,22)=T526;
  g1(35,11)=(-(exp(y(11))/(exp(y(67))+exp(y(11))+exp(y(23)))));
  g1(35,67)=(-(exp(y(67))/(exp(y(67))+exp(y(11))+exp(y(23)))));
  g1(35,73)=1;
  g1(35,23)=(-(exp(y(23))/(exp(y(67))+exp(y(11))+exp(y(23)))));
  g1(36,3)=T459;
  g1(36,46)=T459;
  g1(36,74)=1;
  g1(36,24)=T459;
  g1(37,75)=1;
  g1(37,26)=(-1);
  g1(38,76)=1;
  g1(38,28)=(-1);
  g1(39,77)=1;
  g1(39,30)=(-1);
  g1(40,78)=1;
  g1(40,32)=(-1);
  g1(41,79)=1;
  g1(41,34)=(-1);
  g1(42,80)=1;
  g1(42,36)=(-1);
  g1(43,81)=1;
  g1(43,38)=(-1);
  g1(44,5)=(-1);
  g1(44,82)=1;
  g1(45,6)=(-1);
  g1(45,83)=1;
  g1(46,7)=(-1);
  g1(46,84)=1;
  g1(47,8)=(-1);
  g1(47,85)=1;
  g1(48,11)=(-1);
  g1(48,86)=1;
  g1(49,3)=(-1);
  g1(49,87)=1;
  g1(50,12)=(-1);
  g1(50,88)=1;
  g1(51,25)=(-1);
  g1(51,89)=1;
  g1(52,13)=(-1);
  g1(52,90)=1;
  g1(53,27)=(-1);
  g1(53,91)=1;
  g1(54,14)=(-1);
  g1(54,92)=1;
  g1(55,29)=(-1);
  g1(55,93)=1;
  g1(56,15)=(-1);
  g1(56,94)=1;
  g1(57,31)=(-1);
  g1(57,95)=1;
  g1(58,16)=(-1);
  g1(58,96)=1;
  g1(59,33)=(-1);
  g1(59,97)=1;
  g1(60,17)=(-1);
  g1(60,98)=1;
  g1(61,35)=(-1);
  g1(61,99)=1;
  g1(62,18)=(-1);
  g1(62,100)=1;
  g1(63,37)=(-1);
  g1(63,101)=1;
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],63,12544);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],63,1404928);
end
end
