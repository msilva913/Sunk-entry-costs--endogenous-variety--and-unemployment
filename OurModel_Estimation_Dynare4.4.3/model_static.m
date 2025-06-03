function [residual, g1, g2] = model_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                     columns: variables in declaration order
%                                                     rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 58, 1);

%
% Model equations
%

T47 = 1-exp(y(28));
T50 = y(5)*params(12)*T47/y(5);
T68 = T47*y(5)*params(12)/y(5);
T86 = y(11)/params(19);
T96 = y(1)/T47;
lhs =y(1);
rhs =y(16)/y(15);
residual(1)= lhs-rhs;
lhs =y(2);
rhs =params(16)*y(1)^(1-params(13));
residual(2)= lhs-rhs;
lhs =y(3);
rhs =params(16)*y(1)^(-params(13));
residual(3)= lhs-rhs;
lhs =y(4);
rhs =y(14)^(1/(params(5)-1));
residual(4)= lhs-rhs;
lhs =y(5);
rhs =y(21)^(-params(14));
residual(5)= lhs-rhs;
lhs =y(6);
rhs =(params(5)-1)*y(4)*params(17)/params(5);
residual(6)= lhs-rhs;
lhs =y(6);
rhs =T50*(y(6)+y(7));
residual(7)= lhs-rhs;
lhs =y(8);
rhs =y(18)*exp(y(27))/params(17);
residual(8)= lhs-rhs;
lhs =params(18)+y(9)/y(3);
rhs =T68*(y(12)-y(13)-y(9)+(params(18)+y(9)/y(3))*(1-exp(y(29))));
residual(9)= lhs-rhs;
lhs =y(9);
rhs =y(10)-T68*y(10);
residual(10)= lhs-rhs;
lhs =y(10);
rhs =T86^params(3);
residual(11)= lhs-rhs;
lhs =y(12);
rhs =(params(5)-1)*y(4)*exp(y(27))/params(5);
residual(12)= lhs-rhs;
lhs =y(13);
rhs =params(20)*(y(12)-y(9)+T96*(y(9)+y(3)*params(18)))+(1-params(20))*params(1)*(y(13));
residual(13)= lhs-rhs;
lhs =y(14);
rhs =T47*(y(14)+y(8));
residual(14)= lhs-rhs;
lhs =y(15);
rhs =y(15)*(1-y(2)*T47)+y(26)*(1-y(15));
residual(15)= lhs-rhs;
lhs =y(16);
rhs =y(11)+T47*(y(16)*(1-y(3))+exp(y(29))*(1-y(15)));
residual(16)= lhs-rhs;
lhs =y(17);
rhs =1-y(15);
residual(17)= lhs-rhs;
lhs =y(17);
rhs =y(18)+y(19);
residual(18)= lhs-rhs;
lhs =y(22);
rhs =y(4)*exp(y(27))*y(19);
residual(19)= lhs-rhs;
lhs =y(20);
rhs =params(19)/(1+params(3))*T86^(1+params(3));
residual(20)= lhs-rhs;
lhs =y(22);
rhs =y(21)+y(20)+y(3)*y(16)*params(18);
residual(21)= lhs-rhs;
lhs =y(23);
rhs =y(22)+y(6)*y(8);
residual(22)= lhs-rhs;
lhs =y(23);
rhs =y(12)*y(17)+y(14)*y(7);
residual(23)= lhs-rhs;
lhs =y(24);
rhs =y(13)*y(17)/y(23);
residual(24)= lhs-rhs;
lhs =y(25);
rhs =y(23)/(y(4)*y(17));
residual(25)= lhs-rhs;
lhs =y(26);
rhs =1-T47*(1-exp(y(29)));
residual(26)= lhs-rhs;
lhs =y(27);
rhs =y(27)*params(6)+(1-params(6))*(y(27))-params(9)*x(1);
residual(27)= lhs-rhs;
lhs =y(28);
rhs =y(28)*params(7)+(1-params(7))*(y(28))+params(10)*x(2);
residual(28)= lhs-rhs;
lhs =y(29);
rhs =y(29)*params(8)+(1-params(8))*(y(29))-params(11)*x(3);
residual(29)= lhs-rhs;
lhs =y(30);
rhs =log(0.3333333333333333*(y(15)+y(15)+y(42)));
residual(30)= lhs-rhs;
lhs =y(31);
rhs =log(0.3333333333333333*(y(16)+y(16)+y(43)));
residual(31)= lhs-rhs;
lhs =y(32);
rhs =y(31)-y(30);
residual(32)= lhs-rhs;
lhs =y(33);
rhs =log(0.3333333333333333*(y(25)+y(25)+y(44)));
residual(33)= lhs-rhs;
lhs =y(34);
rhs =log(0.3333333333333333*(y(26)+y(26)+y(45)));
residual(34)= lhs-rhs;
lhs =y(35);
rhs =log(y(8)+y(8)+y(46));
residual(35)= lhs-rhs;
lhs =y(36);
rhs =y(48);
residual(36)= lhs-rhs;
lhs =y(37);
rhs =y(50);
residual(37)= lhs-rhs;
lhs =y(38);
rhs =y(52);
residual(38)= lhs-rhs;
lhs =y(39);
rhs =y(54);
residual(39)= lhs-rhs;
lhs =y(40);
rhs =y(56);
residual(40)= lhs-rhs;
lhs =y(41);
rhs =y(58);
residual(41)= lhs-rhs;
lhs =y(42);
rhs =y(15);
residual(42)= lhs-rhs;
lhs =y(43);
rhs =y(16);
residual(43)= lhs-rhs;
lhs =y(44);
rhs =y(25);
residual(44)= lhs-rhs;
lhs =y(45);
rhs =y(26);
residual(45)= lhs-rhs;
lhs =y(46);
rhs =y(8);
residual(46)= lhs-rhs;
lhs =y(47);
rhs =y(30);
residual(47)= lhs-rhs;
lhs =y(48);
rhs =y(47);
residual(48)= lhs-rhs;
lhs =y(49);
rhs =y(31);
residual(49)= lhs-rhs;
lhs =y(50);
rhs =y(49);
residual(50)= lhs-rhs;
lhs =y(51);
rhs =y(32);
residual(51)= lhs-rhs;
lhs =y(52);
rhs =y(51);
residual(52)= lhs-rhs;
lhs =y(53);
rhs =y(33);
residual(53)= lhs-rhs;
lhs =y(54);
rhs =y(53);
residual(54)= lhs-rhs;
lhs =y(55);
rhs =y(34);
residual(55)= lhs-rhs;
lhs =y(56);
rhs =y(55);
residual(56)= lhs-rhs;
lhs =y(57);
rhs =y(35);
residual(57)= lhs-rhs;
lhs =y(58);
rhs =y(57);
residual(58)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(58, 58);

  %
  % Jacobian matrix
  %

  g1(1,1)=1;
  g1(1,15)=(-((-y(16))/(y(15)*y(15))));
  g1(1,16)=(-(1/y(15)));
  g1(2,1)=(-(params(16)*getPowerDeriv(y(1),1-params(13),1)));
  g1(2,2)=1;
  g1(3,1)=(-(params(16)*getPowerDeriv(y(1),(-params(13)),1)));
  g1(3,3)=1;
  g1(4,4)=1;
  g1(4,14)=(-(getPowerDeriv(y(14),1/(params(5)-1),1)));
  g1(5,5)=1;
  g1(5,21)=(-(getPowerDeriv(y(21),(-params(14)),1)));
  g1(6,4)=(-((params(5)-1)*params(17)/params(5)));
  g1(6,6)=1;
  g1(7,6)=1-T50;
  g1(7,7)=(-T50);
  g1(7,28)=(-((y(6)+y(7))*y(5)*params(12)*(-exp(y(28)))/y(5)));
  g1(8,8)=1;
  g1(8,18)=(-(exp(y(27))/params(17)));
  g1(8,27)=(-(y(18)*exp(y(27))/params(17)));
  g1(9,3)=(-y(9))/(y(3)*y(3))-T68*(1-exp(y(29)))*(-y(9))/(y(3)*y(3));
  g1(9,9)=1/y(3)-T68*((-1)+(1-exp(y(29)))*1/y(3));
  g1(9,12)=(-T68);
  g1(9,13)=T68;
  g1(9,28)=(-((y(12)-y(13)-y(9)+(params(18)+y(9)/y(3))*(1-exp(y(29))))*y(5)*params(12)/y(5)*(-exp(y(28)))));
  g1(9,29)=(-(T68*(params(18)+y(9)/y(3))*(-exp(y(29)))));
  g1(10,9)=1;
  g1(10,10)=(-(1-T68));
  g1(10,28)=y(10)*y(5)*params(12)/y(5)*(-exp(y(28)));
  g1(11,10)=1;
  g1(11,11)=(-(1/params(19)*getPowerDeriv(T86,params(3),1)));
  g1(12,4)=(-((params(5)-1)*exp(y(27))/params(5)));
  g1(12,12)=1;
  g1(12,27)=(-((params(5)-1)*y(4)*exp(y(27))/params(5)));
  g1(13,1)=(-(params(20)*(y(9)+y(3)*params(18))*1/T47));
  g1(13,3)=(-(params(20)*params(18)*T96));
  g1(13,9)=(-(params(20)*((-1)+T96)));
  g1(13,12)=(-params(20));
  g1(13,13)=1-(1-params(20))*params(1);
  g1(13,28)=(-(params(20)*(y(9)+y(3)*params(18))*(-(y(1)*(-exp(y(28)))))/(T47*T47)));
  g1(14,8)=(-T47);
  g1(14,14)=1-T47;
  g1(14,28)=(-((y(14)+y(8))*(-exp(y(28)))));
  g1(15,2)=(-(y(15)*(-T47)));
  g1(15,15)=1-(1-y(2)*T47-y(26));
  g1(15,26)=(-(1-y(15)));
  g1(15,28)=(-(y(15)*(-(y(2)*(-exp(y(28)))))));
  g1(16,3)=(-(T47*(-y(16))));
  g1(16,11)=(-1);
  g1(16,15)=(-(T47*(-exp(y(29)))));
  g1(16,16)=1-T47*(1-y(3));
  g1(16,28)=(-((y(16)*(1-y(3))+exp(y(29))*(1-y(15)))*(-exp(y(28)))));
  g1(16,29)=(-(T47*exp(y(29))*(1-y(15))));
  g1(17,15)=1;
  g1(17,17)=1;
  g1(18,17)=1;
  g1(18,18)=(-1);
  g1(18,19)=(-1);
  g1(19,4)=(-(exp(y(27))*y(19)));
  g1(19,19)=(-(y(4)*exp(y(27))));
  g1(19,22)=1;
  g1(19,27)=(-(y(4)*exp(y(27))*y(19)));
  g1(20,11)=(-(params(19)/(1+params(3))*1/params(19)*getPowerDeriv(T86,1+params(3),1)));
  g1(20,20)=1;
  g1(21,3)=(-(y(16)*params(18)));
  g1(21,16)=(-(y(3)*params(18)));
  g1(21,20)=(-1);
  g1(21,21)=(-1);
  g1(21,22)=1;
  g1(22,6)=(-y(8));
  g1(22,8)=(-y(6));
  g1(22,22)=(-1);
  g1(22,23)=1;
  g1(23,7)=(-y(14));
  g1(23,12)=(-y(17));
  g1(23,14)=(-y(7));
  g1(23,17)=(-y(12));
  g1(23,23)=1;
  g1(24,13)=(-(y(17)/y(23)));
  g1(24,17)=(-(y(13)/y(23)));
  g1(24,23)=(-((-(y(13)*y(17)))/(y(23)*y(23))));
  g1(24,24)=1;
  g1(25,4)=(-((-(y(17)*y(23)))/(y(4)*y(17)*y(4)*y(17))));
  g1(25,17)=(-((-(y(4)*y(23)))/(y(4)*y(17)*y(4)*y(17))));
  g1(25,23)=(-(1/(y(4)*y(17))));
  g1(25,25)=1;
  g1(26,26)=1;
  g1(26,28)=(1-exp(y(29)))*(-exp(y(28)));
  g1(26,29)=T47*(-exp(y(29)));
  g1(27,27)=1-(params(6)+1-params(6));
  g1(28,28)=1-(params(7)+1-params(7));
  g1(29,29)=1-(params(8)+1-params(8));
  g1(30,15)=(-(0.6666666666666666/(0.3333333333333333*(y(15)+y(15)+y(42)))));
  g1(30,30)=1;
  g1(30,42)=(-(0.3333333333333333/(0.3333333333333333*(y(15)+y(15)+y(42)))));
  g1(31,16)=(-(0.6666666666666666/(0.3333333333333333*(y(16)+y(16)+y(43)))));
  g1(31,31)=1;
  g1(31,43)=(-(0.3333333333333333/(0.3333333333333333*(y(16)+y(16)+y(43)))));
  g1(32,30)=1;
  g1(32,31)=(-1);
  g1(32,32)=1;
  g1(33,25)=(-(0.6666666666666666/(0.3333333333333333*(y(25)+y(25)+y(44)))));
  g1(33,33)=1;
  g1(33,44)=(-(0.3333333333333333/(0.3333333333333333*(y(25)+y(25)+y(44)))));
  g1(34,26)=(-(0.6666666666666666/(0.3333333333333333*(y(26)+y(26)+y(45)))));
  g1(34,34)=1;
  g1(34,45)=(-(0.3333333333333333/(0.3333333333333333*(y(26)+y(26)+y(45)))));
  g1(35,8)=(-(2/(y(8)+y(8)+y(46))));
  g1(35,35)=1;
  g1(35,46)=(-(1/(y(8)+y(8)+y(46))));
  g1(36,36)=1;
  g1(36,48)=(-1);
  g1(37,37)=1;
  g1(37,50)=(-1);
  g1(38,38)=1;
  g1(38,52)=(-1);
  g1(39,39)=1;
  g1(39,54)=(-1);
  g1(40,40)=1;
  g1(40,56)=(-1);
  g1(41,41)=1;
  g1(41,58)=(-1);
  g1(42,15)=(-1);
  g1(42,42)=1;
  g1(43,16)=(-1);
  g1(43,43)=1;
  g1(44,25)=(-1);
  g1(44,44)=1;
  g1(45,26)=(-1);
  g1(45,45)=1;
  g1(46,8)=(-1);
  g1(46,46)=1;
  g1(47,30)=(-1);
  g1(47,47)=1;
  g1(48,47)=(-1);
  g1(48,48)=1;
  g1(49,31)=(-1);
  g1(49,49)=1;
  g1(50,49)=(-1);
  g1(50,50)=1;
  g1(51,32)=(-1);
  g1(51,51)=1;
  g1(52,51)=(-1);
  g1(52,52)=1;
  g1(53,33)=(-1);
  g1(53,53)=1;
  g1(54,53)=(-1);
  g1(54,54)=1;
  g1(55,34)=(-1);
  g1(55,55)=1;
  g1(56,55)=(-1);
  g1(56,56)=1;
  g1(57,35)=(-1);
  g1(57,57)=1;
  g1(58,57)=(-1);
  g1(58,58)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],58,3364);
end
end
