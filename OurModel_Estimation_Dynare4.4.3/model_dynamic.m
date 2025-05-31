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

residual = zeros(57, 1);
T47 = 1-exp(y(59));
T70 = T47*params(12)*y(91)/y(37);
T83 = y(96)-y(97)-y(94)+(1-exp(y(60)))*(params(18)+y(94)/y(90));
T93 = y(43)/params(19);
T105 = y(33)/T47;
lhs =y(33);
rhs =y(48)/y(47);
residual(1)= lhs-rhs;
lhs =y(34);
rhs =params(16)*y(33)^(1-params(13));
residual(2)= lhs-rhs;
lhs =y(35);
rhs =params(16)*y(33)^(-params(13));
residual(3)= lhs-rhs;
lhs =y(36);
rhs =y(46)^(1/(params(5)-1));
residual(4)= lhs-rhs;
lhs =y(37);
rhs =y(53)^(-params(14));
residual(5)= lhs-rhs;
lhs =y(38);
rhs =(params(5)-1)*y(36)*params(17)/params(5);
residual(6)= lhs-rhs;
lhs =y(38);
rhs =params(12)*T47*y(91)/y(37)*(y(92)+y(93));
residual(7)= lhs-rhs;
lhs =y(40);
rhs =y(50)*exp(y(58))/params(17);
residual(8)= lhs-rhs;
lhs =params(18)+y(41)/y(35);
rhs =T70*T83;
residual(9)= lhs-rhs;
lhs =y(41);
rhs =y(42)-T70*y(95);
residual(10)= lhs-rhs;
lhs =y(42);
rhs =T93^params(3);
residual(11)= lhs-rhs;
lhs =y(44);
rhs =(params(5)-1)*y(36)*exp(y(58))/params(5);
residual(12)= lhs-rhs;
lhs =y(45);
rhs =params(20)*(y(44)-y(41)+T105*(y(41)+y(35)*params(18)))+(1-params(20))*params(1)*(steady_state(13));
residual(13)= lhs-rhs;
lhs =y(46);
rhs =T47*(y(4)+y(3));
residual(14)= lhs-rhs;
lhs =y(47);
rhs =(1-T47*y(1))*y(5)+y(57)*(1-y(5));
residual(15)= lhs-rhs;
lhs =y(48);
rhs =y(43)+T47*((1-y(2))*y(6)+exp(y(60))*(1-y(5)));
residual(16)= lhs-rhs;
lhs =y(49);
rhs =1-y(47);
residual(17)= lhs-rhs;
lhs =y(49);
rhs =y(50)+y(51);
residual(18)= lhs-rhs;
lhs =y(54);
rhs =y(36)*exp(y(58))*y(51);
residual(19)= lhs-rhs;
lhs =y(52);
rhs =params(19)/(1+params(3))*T93^(1+params(3));
residual(20)= lhs-rhs;
lhs =y(54);
rhs =y(53)+y(52)+y(35)*y(48)*params(18);
residual(21)= lhs-rhs;
lhs =y(55);
rhs =y(54)+y(38)*y(40);
residual(22)= lhs-rhs;
lhs =y(55);
rhs =y(44)*y(49)+y(46)*y(39);
residual(23)= lhs-rhs;
lhs =y(56);
rhs =y(45)*y(49)/y(55);
residual(24)= lhs-rhs;
lhs =y(57);
rhs =1-T47*(1-exp(y(60)));
residual(25)= lhs-rhs;
lhs =y(58);
rhs =params(6)*y(7)+(1-params(6))*(steady_state(26))-params(9)*x(it_, 1);
residual(26)= lhs-rhs;
lhs =y(59);
rhs =params(7)*y(8)+(1-params(7))*(steady_state(27))+params(10)*x(it_, 2);
residual(27)= lhs-rhs;
lhs =y(60);
rhs =params(8)*y(9)+(1-params(8))*(steady_state(28))-params(11)*x(it_, 3);
residual(28)= lhs-rhs;
lhs =y(61);
rhs =0.3333333333333333*(y(47)+y(5)+y(16));
residual(29)= lhs-rhs;
lhs =y(62);
rhs =0.3333333333333333*(y(48)+y(6)+y(17));
residual(30)= lhs-rhs;
lhs =y(63);
rhs =y(62)/y(61);
residual(31)= lhs-rhs;
lhs =y(64);
rhs =0.3333333333333333*(exp(y(58))+exp(y(7))+exp(y(18)));
residual(32)= lhs-rhs;
lhs =y(65);
rhs =0.3333333333333333*(exp(y(60))+exp(y(9))+exp(y(19)));
residual(33)= lhs-rhs;
lhs =y(66);
rhs =(y(40)+y(3)+y(20))/(steady_state(14));
residual(34)= lhs-rhs;
lhs =y(67);
rhs =y(22);
residual(35)= lhs-rhs;
lhs =y(68);
rhs =y(24);
residual(36)= lhs-rhs;
lhs =y(69);
rhs =y(26);
residual(37)= lhs-rhs;
lhs =y(70);
rhs =y(28);
residual(38)= lhs-rhs;
lhs =y(71);
rhs =y(30);
residual(39)= lhs-rhs;
lhs =y(72);
rhs =y(32);
residual(40)= lhs-rhs;
lhs =y(73);
rhs =y(5);
residual(41)= lhs-rhs;
lhs =y(74);
rhs =y(6);
residual(42)= lhs-rhs;
lhs =y(75);
rhs =y(7);
residual(43)= lhs-rhs;
lhs =y(76);
rhs =y(9);
residual(44)= lhs-rhs;
lhs =y(77);
rhs =y(3);
residual(45)= lhs-rhs;
lhs =y(78);
rhs =y(10);
residual(46)= lhs-rhs;
lhs =y(79);
rhs =y(21);
residual(47)= lhs-rhs;
lhs =y(80);
rhs =y(11);
residual(48)= lhs-rhs;
lhs =y(81);
rhs =y(23);
residual(49)= lhs-rhs;
lhs =y(82);
rhs =y(12);
residual(50)= lhs-rhs;
lhs =y(83);
rhs =y(25);
residual(51)= lhs-rhs;
lhs =y(84);
rhs =y(13);
residual(52)= lhs-rhs;
lhs =y(85);
rhs =y(27);
residual(53)= lhs-rhs;
lhs =y(86);
rhs =y(14);
residual(54)= lhs-rhs;
lhs =y(87);
rhs =y(29);
residual(55)= lhs-rhs;
lhs =y(88);
rhs =y(15);
residual(56)= lhs-rhs;
lhs =y(89);
rhs =y(31);
residual(57)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(57, 100);

  %
  % Jacobian matrix
  %

  g1(1,33)=1;
  g1(1,47)=(-((-y(48))/(y(47)*y(47))));
  g1(1,48)=(-(1/y(47)));
  g1(2,33)=(-(params(16)*getPowerDeriv(y(33),1-params(13),1)));
  g1(2,34)=1;
  g1(3,33)=(-(params(16)*getPowerDeriv(y(33),(-params(13)),1)));
  g1(3,35)=1;
  g1(4,36)=1;
  g1(4,46)=(-(getPowerDeriv(y(46),1/(params(5)-1),1)));
  g1(5,37)=1;
  g1(5,53)=(-(getPowerDeriv(y(53),(-params(14)),1)));
  g1(6,36)=(-((params(5)-1)*params(17)/params(5)));
  g1(6,38)=1;
  g1(7,37)=(-((y(92)+y(93))*(-(params(12)*T47*y(91)))/(y(37)*y(37))));
  g1(7,91)=(-((y(92)+y(93))*params(12)*T47/y(37)));
  g1(7,38)=1;
  g1(7,92)=(-(params(12)*T47*y(91)/y(37)));
  g1(7,93)=(-(params(12)*T47*y(91)/y(37)));
  g1(7,59)=(-((y(92)+y(93))*y(91)*params(12)*(-exp(y(59)))/y(37)));
  g1(8,40)=1;
  g1(8,50)=(-(exp(y(58))/params(17)));
  g1(8,58)=(-(y(50)*exp(y(58))/params(17)));
  g1(9,35)=(-y(41))/(y(35)*y(35));
  g1(9,90)=(-(T70*(1-exp(y(60)))*(-y(94))/(y(90)*y(90))));
  g1(9,37)=(-(T83*T47*(-(params(12)*y(91)))/(y(37)*y(37))));
  g1(9,91)=(-(T83*T47*params(12)/y(37)));
  g1(9,41)=1/y(35);
  g1(9,94)=(-(T70*((-1)+(1-exp(y(60)))*1/y(90))));
  g1(9,96)=(-T70);
  g1(9,97)=T70;
  g1(9,59)=(-(T83*params(12)*y(91)/y(37)*(-exp(y(59)))));
  g1(9,60)=(-(T70*(params(18)+y(94)/y(90))*(-exp(y(60)))));
  g1(10,37)=y(95)*T47*(-(params(12)*y(91)))/(y(37)*y(37));
  g1(10,91)=y(95)*T47*params(12)/y(37);
  g1(10,41)=1;
  g1(10,42)=(-1);
  g1(10,95)=T70;
  g1(10,59)=y(95)*params(12)*y(91)/y(37)*(-exp(y(59)));
  g1(11,42)=1;
  g1(11,43)=(-(1/params(19)*getPowerDeriv(T93,params(3),1)));
  g1(12,36)=(-((params(5)-1)*exp(y(58))/params(5)));
  g1(12,44)=1;
  g1(12,58)=(-((params(5)-1)*y(36)*exp(y(58))/params(5)));
  g1(13,33)=(-(params(20)*(y(41)+y(35)*params(18))*1/T47));
  g1(13,35)=(-(params(20)*params(18)*T105));
  g1(13,41)=(-(params(20)*((-1)+T105)));
  g1(13,44)=(-params(20));
  g1(13,45)=1;
  g1(13,59)=(-(params(20)*(y(41)+y(35)*params(18))*(-(y(33)*(-exp(y(59)))))/(T47*T47)));
  g1(14,3)=(-T47);
  g1(14,4)=(-T47);
  g1(14,46)=1;
  g1(14,59)=(-((y(4)+y(3))*(-exp(y(59)))));
  g1(15,1)=(-(y(5)*(-T47)));
  g1(15,5)=(-(1-T47*y(1)-y(57)));
  g1(15,47)=1;
  g1(15,57)=(-(1-y(5)));
  g1(15,59)=(-(y(5)*(-(y(1)*(-exp(y(59)))))));
  g1(16,2)=(-(T47*(-y(6))));
  g1(16,43)=(-1);
  g1(16,5)=(-(T47*(-exp(y(60)))));
  g1(16,6)=(-(T47*(1-y(2))));
  g1(16,48)=1;
  g1(16,59)=(-(((1-y(2))*y(6)+exp(y(60))*(1-y(5)))*(-exp(y(59)))));
  g1(16,60)=(-(T47*exp(y(60))*(1-y(5))));
  g1(17,47)=1;
  g1(17,49)=1;
  g1(18,49)=1;
  g1(18,50)=(-1);
  g1(18,51)=(-1);
  g1(19,36)=(-(exp(y(58))*y(51)));
  g1(19,51)=(-(y(36)*exp(y(58))));
  g1(19,54)=1;
  g1(19,58)=(-(y(36)*exp(y(58))*y(51)));
  g1(20,43)=(-(params(19)/(1+params(3))*1/params(19)*getPowerDeriv(T93,1+params(3),1)));
  g1(20,52)=1;
  g1(21,35)=(-(y(48)*params(18)));
  g1(21,48)=(-(y(35)*params(18)));
  g1(21,52)=(-1);
  g1(21,53)=(-1);
  g1(21,54)=1;
  g1(22,38)=(-y(40));
  g1(22,40)=(-y(38));
  g1(22,54)=(-1);
  g1(22,55)=1;
  g1(23,39)=(-y(46));
  g1(23,44)=(-y(49));
  g1(23,46)=(-y(39));
  g1(23,49)=(-y(44));
  g1(23,55)=1;
  g1(24,45)=(-(y(49)/y(55)));
  g1(24,49)=(-(y(45)/y(55)));
  g1(24,55)=(-((-(y(45)*y(49)))/(y(55)*y(55))));
  g1(24,56)=1;
  g1(25,57)=1;
  g1(25,59)=(1-exp(y(60)))*(-exp(y(59)));
  g1(25,60)=T47*(-exp(y(60)));
  g1(26,7)=(-params(6));
  g1(26,58)=1;
  g1(26,98)=params(9);
  g1(27,8)=(-params(7));
  g1(27,59)=1;
  g1(27,99)=(-params(10));
  g1(28,9)=(-params(8));
  g1(28,60)=1;
  g1(28,100)=params(11);
  g1(29,5)=(-0.3333333333333333);
  g1(29,47)=(-0.3333333333333333);
  g1(29,61)=1;
  g1(29,16)=(-0.3333333333333333);
  g1(30,6)=(-0.3333333333333333);
  g1(30,48)=(-0.3333333333333333);
  g1(30,62)=1;
  g1(30,17)=(-0.3333333333333333);
  g1(31,61)=(-((-y(62))/(y(61)*y(61))));
  g1(31,62)=(-(1/y(61)));
  g1(31,63)=1;
  g1(32,7)=(-(0.3333333333333333*exp(y(7))));
  g1(32,58)=(-(exp(y(58))*0.3333333333333333));
  g1(32,64)=1;
  g1(32,18)=(-(0.3333333333333333*exp(y(18))));
  g1(33,9)=(-(0.3333333333333333*exp(y(9))));
  g1(33,60)=(-(exp(y(60))*0.3333333333333333));
  g1(33,65)=1;
  g1(33,19)=(-(0.3333333333333333*exp(y(19))));
  g1(34,3)=(-(1/(steady_state(14))));
  g1(34,40)=(-(1/(steady_state(14))));
  g1(34,66)=1;
  g1(34,20)=(-(1/(steady_state(14))));
  g1(35,67)=1;
  g1(35,22)=(-1);
  g1(36,68)=1;
  g1(36,24)=(-1);
  g1(37,69)=1;
  g1(37,26)=(-1);
  g1(38,70)=1;
  g1(38,28)=(-1);
  g1(39,71)=1;
  g1(39,30)=(-1);
  g1(40,72)=1;
  g1(40,32)=(-1);
  g1(41,5)=(-1);
  g1(41,73)=1;
  g1(42,6)=(-1);
  g1(42,74)=1;
  g1(43,7)=(-1);
  g1(43,75)=1;
  g1(44,9)=(-1);
  g1(44,76)=1;
  g1(45,3)=(-1);
  g1(45,77)=1;
  g1(46,10)=(-1);
  g1(46,78)=1;
  g1(47,21)=(-1);
  g1(47,79)=1;
  g1(48,11)=(-1);
  g1(48,80)=1;
  g1(49,23)=(-1);
  g1(49,81)=1;
  g1(50,12)=(-1);
  g1(50,82)=1;
  g1(51,25)=(-1);
  g1(51,83)=1;
  g1(52,13)=(-1);
  g1(52,84)=1;
  g1(53,27)=(-1);
  g1(53,85)=1;
  g1(54,14)=(-1);
  g1(54,86)=1;
  g1(55,29)=(-1);
  g1(55,87)=1;
  g1(56,15)=(-1);
  g1(56,88)=1;
  g1(57,31)=(-1);
  g1(57,89)=1;
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],57,10000);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],57,1000000);
end
end
