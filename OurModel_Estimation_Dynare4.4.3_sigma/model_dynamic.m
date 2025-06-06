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

residual = zeros(58, 1);
T47 = 1-exp(y(62));
T70 = T47*params(13)*y(94)/y(39);
T83 = y(99)-y(100)-y(97)+(1-exp(y(63)))*(params(18)+y(97)/y(93));
T93 = y(45)/params(19);
T105 = y(35)/T47;
T434 = (-(1/(y(42)+y(3)+y(22))));
T468 = (-(0.3333333333333333/(0.3333333333333333*(y(49)+y(5)+y(18)))));
T476 = (-(0.3333333333333333/(0.3333333333333333*(y(50)+y(6)+y(19)))));
T499 = (-(0.3333333333333333/(0.3333333333333333*(y(59)+y(7)+y(20)))));
T501 = (-(0.3333333333333333/(0.3333333333333333*(y(60)+y(8)+y(21)))));
lhs =y(35);
rhs =y(50)/y(49);
residual(1)= lhs-rhs;
lhs =y(36);
rhs =params(16)*y(35)^(1-params(14));
residual(2)= lhs-rhs;
lhs =y(37);
rhs =params(16)*y(35)^(-params(14));
residual(3)= lhs-rhs;
lhs =y(38);
rhs =y(48)^(1/(params(6)-1));
residual(4)= lhs-rhs;
lhs =y(39);
rhs =y(55)^(-params(1));
residual(5)= lhs-rhs;
lhs =y(40);
rhs =(params(6)-1)*y(38)*params(17)/params(6);
residual(6)= lhs-rhs;
lhs =y(40);
rhs =params(13)*T47*y(94)/y(39)*(y(95)+y(96));
residual(7)= lhs-rhs;
lhs =y(42);
rhs =y(52)*exp(y(61))/params(17);
residual(8)= lhs-rhs;
lhs =params(18)+y(43)/y(37);
rhs =T70*T83;
residual(9)= lhs-rhs;
lhs =y(43);
rhs =y(44)-T70*y(98);
residual(10)= lhs-rhs;
lhs =y(44);
rhs =T93^params(4);
residual(11)= lhs-rhs;
lhs =y(46);
rhs =(params(6)-1)*y(38)*exp(y(61))/params(6);
residual(12)= lhs-rhs;
lhs =y(47);
rhs =params(20)*(y(46)-y(43)+T105*(y(43)+y(37)*params(18)))+(1-params(20))*params(2)*(steady_state(13));
residual(13)= lhs-rhs;
lhs =y(48);
rhs =T47*(y(4)+y(3));
residual(14)= lhs-rhs;
lhs =y(49);
rhs =(1-T47*y(1))*y(5)+y(60)*(1-y(5));
residual(15)= lhs-rhs;
lhs =y(50);
rhs =y(45)+T47*((1-y(2))*y(6)+exp(y(63))*(1-y(5)));
residual(16)= lhs-rhs;
lhs =y(51);
rhs =1-y(49);
residual(17)= lhs-rhs;
lhs =y(51);
rhs =y(52)+y(53);
residual(18)= lhs-rhs;
lhs =y(56);
rhs =y(38)*exp(y(61))*y(53);
residual(19)= lhs-rhs;
lhs =y(54);
rhs =params(19)/(1+params(4))*T93^(1+params(4));
residual(20)= lhs-rhs;
lhs =y(56);
rhs =y(55)+y(54)+y(37)*y(50)*params(18);
residual(21)= lhs-rhs;
lhs =y(57);
rhs =y(56)+y(40)*y(42);
residual(22)= lhs-rhs;
lhs =y(57);
rhs =y(46)*y(51)+y(48)*y(41);
residual(23)= lhs-rhs;
lhs =y(58);
rhs =y(47)*y(51)/y(57);
residual(24)= lhs-rhs;
lhs =y(59);
rhs =y(57)/(y(38)*y(51));
residual(25)= lhs-rhs;
lhs =y(60);
rhs =1-T47*(1-exp(y(63)));
residual(26)= lhs-rhs;
lhs =y(61);
rhs =params(7)*y(9)+(1-params(7))*(steady_state(27))-params(10)*x(it_, 1);
residual(27)= lhs-rhs;
lhs =y(62);
rhs =params(8)*y(10)+(1-params(8))*(steady_state(28))+params(11)*x(it_, 2);
residual(28)= lhs-rhs;
lhs =y(63);
rhs =params(9)*y(11)+(1-params(9))*(steady_state(29))-params(12)*x(it_, 3);
residual(29)= lhs-rhs;
lhs =y(64);
rhs =log(0.3333333333333333*(y(49)+y(5)+y(18)));
residual(30)= lhs-rhs;
lhs =y(65);
rhs =log(0.3333333333333333*(y(50)+y(6)+y(19)));
residual(31)= lhs-rhs;
lhs =y(66);
rhs =y(65)-y(64);
residual(32)= lhs-rhs;
lhs =y(67);
rhs =log(0.3333333333333333*(y(59)+y(7)+y(20)));
residual(33)= lhs-rhs;
lhs =y(68);
rhs =log(0.3333333333333333*(y(60)+y(8)+y(21)));
residual(34)= lhs-rhs;
lhs =y(69);
rhs =log(y(42)+y(3)+y(22));
residual(35)= lhs-rhs;
lhs =y(70);
rhs =y(24);
residual(36)= lhs-rhs;
lhs =y(71);
rhs =y(26);
residual(37)= lhs-rhs;
lhs =y(72);
rhs =y(28);
residual(38)= lhs-rhs;
lhs =y(73);
rhs =y(30);
residual(39)= lhs-rhs;
lhs =y(74);
rhs =y(32);
residual(40)= lhs-rhs;
lhs =y(75);
rhs =y(34);
residual(41)= lhs-rhs;
lhs =y(76);
rhs =y(5);
residual(42)= lhs-rhs;
lhs =y(77);
rhs =y(6);
residual(43)= lhs-rhs;
lhs =y(78);
rhs =y(7);
residual(44)= lhs-rhs;
lhs =y(79);
rhs =y(8);
residual(45)= lhs-rhs;
lhs =y(80);
rhs =y(3);
residual(46)= lhs-rhs;
lhs =y(81);
rhs =y(12);
residual(47)= lhs-rhs;
lhs =y(82);
rhs =y(23);
residual(48)= lhs-rhs;
lhs =y(83);
rhs =y(13);
residual(49)= lhs-rhs;
lhs =y(84);
rhs =y(25);
residual(50)= lhs-rhs;
lhs =y(85);
rhs =y(14);
residual(51)= lhs-rhs;
lhs =y(86);
rhs =y(27);
residual(52)= lhs-rhs;
lhs =y(87);
rhs =y(15);
residual(53)= lhs-rhs;
lhs =y(88);
rhs =y(29);
residual(54)= lhs-rhs;
lhs =y(89);
rhs =y(16);
residual(55)= lhs-rhs;
lhs =y(90);
rhs =y(31);
residual(56)= lhs-rhs;
lhs =y(91);
rhs =y(17);
residual(57)= lhs-rhs;
lhs =y(92);
rhs =y(33);
residual(58)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(58, 103);

  %
  % Jacobian matrix
  %

  g1(1,35)=1;
  g1(1,49)=(-((-y(50))/(y(49)*y(49))));
  g1(1,50)=(-(1/y(49)));
  g1(2,35)=(-(params(16)*getPowerDeriv(y(35),1-params(14),1)));
  g1(2,36)=1;
  g1(3,35)=(-(params(16)*getPowerDeriv(y(35),(-params(14)),1)));
  g1(3,37)=1;
  g1(4,38)=1;
  g1(4,48)=(-(getPowerDeriv(y(48),1/(params(6)-1),1)));
  g1(5,39)=1;
  g1(5,55)=(-(getPowerDeriv(y(55),(-params(1)),1)));
  g1(6,38)=(-((params(6)-1)*params(17)/params(6)));
  g1(6,40)=1;
  g1(7,39)=(-((y(95)+y(96))*(-(params(13)*T47*y(94)))/(y(39)*y(39))));
  g1(7,94)=(-((y(95)+y(96))*params(13)*T47/y(39)));
  g1(7,40)=1;
  g1(7,95)=(-(params(13)*T47*y(94)/y(39)));
  g1(7,96)=(-(params(13)*T47*y(94)/y(39)));
  g1(7,62)=(-((y(95)+y(96))*y(94)*params(13)*(-exp(y(62)))/y(39)));
  g1(8,42)=1;
  g1(8,52)=(-(exp(y(61))/params(17)));
  g1(8,61)=(-(y(52)*exp(y(61))/params(17)));
  g1(9,37)=(-y(43))/(y(37)*y(37));
  g1(9,93)=(-(T70*(1-exp(y(63)))*(-y(97))/(y(93)*y(93))));
  g1(9,39)=(-(T83*T47*(-(params(13)*y(94)))/(y(39)*y(39))));
  g1(9,94)=(-(T83*T47*params(13)/y(39)));
  g1(9,43)=1/y(37);
  g1(9,97)=(-(T70*((-1)+(1-exp(y(63)))*1/y(93))));
  g1(9,99)=(-T70);
  g1(9,100)=T70;
  g1(9,62)=(-(T83*params(13)*y(94)/y(39)*(-exp(y(62)))));
  g1(9,63)=(-(T70*(params(18)+y(97)/y(93))*(-exp(y(63)))));
  g1(10,39)=y(98)*T47*(-(params(13)*y(94)))/(y(39)*y(39));
  g1(10,94)=y(98)*T47*params(13)/y(39);
  g1(10,43)=1;
  g1(10,44)=(-1);
  g1(10,98)=T70;
  g1(10,62)=y(98)*params(13)*y(94)/y(39)*(-exp(y(62)));
  g1(11,44)=1;
  g1(11,45)=(-(1/params(19)*getPowerDeriv(T93,params(4),1)));
  g1(12,38)=(-((params(6)-1)*exp(y(61))/params(6)));
  g1(12,46)=1;
  g1(12,61)=(-((params(6)-1)*y(38)*exp(y(61))/params(6)));
  g1(13,35)=(-(params(20)*(y(43)+y(37)*params(18))*1/T47));
  g1(13,37)=(-(params(20)*params(18)*T105));
  g1(13,43)=(-(params(20)*((-1)+T105)));
  g1(13,46)=(-params(20));
  g1(13,47)=1;
  g1(13,62)=(-(params(20)*(y(43)+y(37)*params(18))*(-(y(35)*(-exp(y(62)))))/(T47*T47)));
  g1(14,3)=(-T47);
  g1(14,4)=(-T47);
  g1(14,48)=1;
  g1(14,62)=(-((y(4)+y(3))*(-exp(y(62)))));
  g1(15,1)=(-(y(5)*(-T47)));
  g1(15,5)=(-(1-T47*y(1)-y(60)));
  g1(15,49)=1;
  g1(15,60)=(-(1-y(5)));
  g1(15,62)=(-(y(5)*(-(y(1)*(-exp(y(62)))))));
  g1(16,2)=(-(T47*(-y(6))));
  g1(16,45)=(-1);
  g1(16,5)=(-(T47*(-exp(y(63)))));
  g1(16,6)=(-(T47*(1-y(2))));
  g1(16,50)=1;
  g1(16,62)=(-(((1-y(2))*y(6)+exp(y(63))*(1-y(5)))*(-exp(y(62)))));
  g1(16,63)=(-(T47*exp(y(63))*(1-y(5))));
  g1(17,49)=1;
  g1(17,51)=1;
  g1(18,51)=1;
  g1(18,52)=(-1);
  g1(18,53)=(-1);
  g1(19,38)=(-(exp(y(61))*y(53)));
  g1(19,53)=(-(y(38)*exp(y(61))));
  g1(19,56)=1;
  g1(19,61)=(-(y(38)*exp(y(61))*y(53)));
  g1(20,45)=(-(params(19)/(1+params(4))*1/params(19)*getPowerDeriv(T93,1+params(4),1)));
  g1(20,54)=1;
  g1(21,37)=(-(y(50)*params(18)));
  g1(21,50)=(-(y(37)*params(18)));
  g1(21,54)=(-1);
  g1(21,55)=(-1);
  g1(21,56)=1;
  g1(22,40)=(-y(42));
  g1(22,42)=(-y(40));
  g1(22,56)=(-1);
  g1(22,57)=1;
  g1(23,41)=(-y(48));
  g1(23,46)=(-y(51));
  g1(23,48)=(-y(41));
  g1(23,51)=(-y(46));
  g1(23,57)=1;
  g1(24,47)=(-(y(51)/y(57)));
  g1(24,51)=(-(y(47)/y(57)));
  g1(24,57)=(-((-(y(47)*y(51)))/(y(57)*y(57))));
  g1(24,58)=1;
  g1(25,38)=(-((-(y(51)*y(57)))/(y(38)*y(51)*y(38)*y(51))));
  g1(25,51)=(-((-(y(38)*y(57)))/(y(38)*y(51)*y(38)*y(51))));
  g1(25,57)=(-(1/(y(38)*y(51))));
  g1(25,59)=1;
  g1(26,60)=1;
  g1(26,62)=(1-exp(y(63)))*(-exp(y(62)));
  g1(26,63)=T47*(-exp(y(63)));
  g1(27,9)=(-params(7));
  g1(27,61)=1;
  g1(27,101)=params(10);
  g1(28,10)=(-params(8));
  g1(28,62)=1;
  g1(28,102)=(-params(11));
  g1(29,11)=(-params(9));
  g1(29,63)=1;
  g1(29,103)=params(12);
  g1(30,5)=T468;
  g1(30,49)=T468;
  g1(30,64)=1;
  g1(30,18)=T468;
  g1(31,6)=T476;
  g1(31,50)=T476;
  g1(31,65)=1;
  g1(31,19)=T476;
  g1(32,64)=1;
  g1(32,65)=(-1);
  g1(32,66)=1;
  g1(33,7)=T499;
  g1(33,59)=T499;
  g1(33,67)=1;
  g1(33,20)=T499;
  g1(34,8)=T501;
  g1(34,60)=T501;
  g1(34,68)=1;
  g1(34,21)=T501;
  g1(35,3)=T434;
  g1(35,42)=T434;
  g1(35,69)=1;
  g1(35,22)=T434;
  g1(36,70)=1;
  g1(36,24)=(-1);
  g1(37,71)=1;
  g1(37,26)=(-1);
  g1(38,72)=1;
  g1(38,28)=(-1);
  g1(39,73)=1;
  g1(39,30)=(-1);
  g1(40,74)=1;
  g1(40,32)=(-1);
  g1(41,75)=1;
  g1(41,34)=(-1);
  g1(42,5)=(-1);
  g1(42,76)=1;
  g1(43,6)=(-1);
  g1(43,77)=1;
  g1(44,7)=(-1);
  g1(44,78)=1;
  g1(45,8)=(-1);
  g1(45,79)=1;
  g1(46,3)=(-1);
  g1(46,80)=1;
  g1(47,12)=(-1);
  g1(47,81)=1;
  g1(48,23)=(-1);
  g1(48,82)=1;
  g1(49,13)=(-1);
  g1(49,83)=1;
  g1(50,25)=(-1);
  g1(50,84)=1;
  g1(51,14)=(-1);
  g1(51,85)=1;
  g1(52,27)=(-1);
  g1(52,86)=1;
  g1(53,15)=(-1);
  g1(53,87)=1;
  g1(54,29)=(-1);
  g1(54,88)=1;
  g1(55,16)=(-1);
  g1(55,89)=1;
  g1(56,31)=(-1);
  g1(56,90)=1;
  g1(57,17)=(-1);
  g1(57,91)=1;
  g1(58,33)=(-1);
  g1(58,92)=1;
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],58,10609);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],58,1092727);
end
end
