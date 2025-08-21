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

residual = zeros( 65, 1);

%
% Model equations
%

mu__ = params(13)/(params(13)-1);
w_ss__ = (y(13));
T47 = y(5)*params(16)*(1-y(30))/y(5);
T66 = params(17)+y(9)/y(3);
T228 = 0.3333333333333333*(y(2)*(1-y(30))+y(2)*(1-y(30))+(1-y(48))*y(49));
lhs =y(1);
rhs =y(16)/y(15);
residual(1)= lhs-rhs;
lhs =y(2);
rhs =params(19)*y(1)^(1-params(11));
residual(2)= lhs-rhs;
lhs =y(3);
rhs =params(19)*y(1)^(-params(11));
residual(3)= lhs-rhs;
lhs =y(4);
rhs =y(14)^params(14);
residual(4)= lhs-rhs;
lhs =y(5);
rhs =y(22)^(-params(1));
residual(5)= lhs-rhs;
lhs =y(6);
rhs =T47*(y(6)+y(7));
residual(6)= lhs-rhs;
lhs =y(6)*y(8);
rhs =y(12)*y(19);
residual(7)= lhs-rhs;
lhs =y(8)*params(18);
rhs =y(19)*y(28);
residual(8)= lhs-rhs;
lhs =T66;
rhs =T47*(y(12)-y(13)-y(9)+T66*(1-y(29)));
residual(9)= lhs-rhs;
lhs =y(9);
rhs =y(10)-T47*y(10);
residual(10)= lhs-rhs;
lhs =y(10);
rhs =params(20)*(y(11)/params(12))^params(3);
residual(11)= lhs-rhs;
lhs =y(12);
rhs =y(4)*y(28)/mu__;
residual(12)= lhs-rhs;
lhs =y(13);
rhs =params(21)*(y(12)-y(9)+y(2)*T66)+(1-params(21))*params(2)*w_ss__;
residual(13)= lhs-rhs;
lhs =y(14);
rhs =(1-y(30))*(y(14)+y(8));
residual(14)= lhs-rhs;
lhs =y(15);
rhs =y(15)*(1-y(2)*(1-y(30)))+y(17)*(1-y(15));
residual(15)= lhs-rhs;
lhs =y(16);
rhs =y(11)+(1-y(30))*(y(16)*(1-y(3))+y(29)*(1-y(15)));
residual(16)= lhs-rhs;
lhs =y(17);
rhs =1-(1-y(30))*(1-y(29));
residual(17)= lhs-rhs;
lhs =y(18);
rhs =1-y(15);
residual(18)= lhs-rhs;
lhs =y(18);
rhs =y(19)+y(20);
residual(19)= lhs-rhs;
lhs =y(23);
rhs =y(4)*y(28)*y(20);
residual(20)= lhs-rhs;
lhs =y(21);
rhs =y(10)*y(11)/(1+params(3))+y(3)*y(16)*params(17);
residual(21)= lhs-rhs;
lhs =y(23);
rhs =y(22)+y(21);
residual(22)= lhs-rhs;
lhs =y(24);
rhs =y(6)*y(8)+y(23);
residual(23)= lhs-rhs;
lhs =y(24);
rhs =y(12)*y(18)+y(14)*y(7);
residual(24)= lhs-rhs;
lhs =y(25);
rhs =y(24)-y(21);
residual(25)= lhs-rhs;
lhs =y(26);
rhs =y(13)*y(18)/y(25);
residual(26)= lhs-rhs;
lhs =y(27);
rhs =y(25)/(y(4)*y(18));
residual(27)= lhs-rhs;
lhs =log(y(28));
rhs =log(y(28))*params(5)+(1-params(5))*log((y(28)))-params(8)*x(1);
residual(28)= lhs-rhs;
lhs =log(y(29));
rhs =log(y(29))*params(6)+(1-params(6))*log((y(29)))+params(9)*x(2);
residual(29)= lhs-rhs;
lhs =log(y(30));
rhs =log(y(30))*params(7)+(1-params(7))*log((y(30)))+params(10)*x(3);
residual(30)= lhs-rhs;
lhs =y(31);
rhs =log(0.3333333333333333*(y(15)+y(15)+y(45)));
residual(31)= lhs-rhs;
lhs =y(32);
rhs =log(0.3333333333333333*(y(16)+y(16)+y(46)));
residual(32)= lhs-rhs;
lhs =y(33);
rhs =log(0.3333333333333333*(y(17)+y(17)+y(47)));
residual(33)= lhs-rhs;
lhs =y(34);
rhs =log(T228);
residual(34)= lhs-rhs;
lhs =y(35);
rhs =log(y(48)+y(30)+y(30));
residual(35)= lhs-rhs;
lhs =y(36);
rhs =log(0.3333333333333333*(y(8)+y(8)+y(50)));
residual(36)= lhs-rhs;
lhs =y(37);
rhs =log(0.3333333333333333*(y(27)+y(27)+y(51)));
residual(37)= lhs-rhs;
lhs =y(38);
rhs =y(53);
residual(38)= lhs-rhs;
lhs =y(39);
rhs =y(55);
residual(39)= lhs-rhs;
lhs =y(40);
rhs =y(57);
residual(40)= lhs-rhs;
lhs =y(41);
rhs =y(59);
residual(41)= lhs-rhs;
lhs =y(42);
rhs =y(61);
residual(42)= lhs-rhs;
lhs =y(43);
rhs =y(63);
residual(43)= lhs-rhs;
lhs =y(44);
rhs =y(65);
residual(44)= lhs-rhs;
lhs =y(45);
rhs =y(15);
residual(45)= lhs-rhs;
lhs =y(46);
rhs =y(16);
residual(46)= lhs-rhs;
lhs =y(47);
rhs =y(17);
residual(47)= lhs-rhs;
lhs =y(48);
rhs =y(30);
residual(48)= lhs-rhs;
lhs =y(49);
rhs =y(2);
residual(49)= lhs-rhs;
lhs =y(50);
rhs =y(8);
residual(50)= lhs-rhs;
lhs =y(51);
rhs =y(27);
residual(51)= lhs-rhs;
lhs =y(52);
rhs =y(31);
residual(52)= lhs-rhs;
lhs =y(53);
rhs =y(52);
residual(53)= lhs-rhs;
lhs =y(54);
rhs =y(32);
residual(54)= lhs-rhs;
lhs =y(55);
rhs =y(54);
residual(55)= lhs-rhs;
lhs =y(56);
rhs =y(33);
residual(56)= lhs-rhs;
lhs =y(57);
rhs =y(56);
residual(57)= lhs-rhs;
lhs =y(58);
rhs =y(34);
residual(58)= lhs-rhs;
lhs =y(59);
rhs =y(58);
residual(59)= lhs-rhs;
lhs =y(60);
rhs =y(35);
residual(60)= lhs-rhs;
lhs =y(61);
rhs =y(60);
residual(61)= lhs-rhs;
lhs =y(62);
rhs =y(36);
residual(62)= lhs-rhs;
lhs =y(63);
rhs =y(62);
residual(63)= lhs-rhs;
lhs =y(64);
rhs =y(37);
residual(64)= lhs-rhs;
lhs =y(65);
rhs =y(64);
residual(65)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(65, 65);

  %
  % Jacobian matrix
  %

  g1(1,1)=1;
  g1(1,15)=(-((-y(16))/(y(15)*y(15))));
  g1(1,16)=(-(1/y(15)));
  g1(2,1)=(-(params(19)*getPowerDeriv(y(1),1-params(11),1)));
  g1(2,2)=1;
  g1(3,1)=(-(params(19)*getPowerDeriv(y(1),(-params(11)),1)));
  g1(3,3)=1;
  g1(4,4)=1;
  g1(4,14)=(-(getPowerDeriv(y(14),params(14),1)));
  g1(5,5)=1;
  g1(5,22)=(-(getPowerDeriv(y(22),(-params(1)),1)));
  g1(6,6)=1-T47;
  g1(6,7)=(-T47);
  g1(6,30)=(-((y(6)+y(7))*y(5)*(-params(16))/y(5)));
  g1(7,6)=y(8);
  g1(7,8)=y(6);
  g1(7,12)=(-y(19));
  g1(7,19)=(-y(12));
  g1(8,8)=params(18);
  g1(8,19)=(-y(28));
  g1(8,28)=(-y(19));
  g1(9,3)=(-y(9))/(y(3)*y(3))-T47*(1-y(29))*(-y(9))/(y(3)*y(3));
  g1(9,9)=1/y(3)-T47*((-1)+(1-y(29))*1/y(3));
  g1(9,12)=(-T47);
  g1(9,13)=T47;
  g1(9,29)=(-(T47*(-T66)));
  g1(9,30)=(-((y(12)-y(13)-y(9)+T66*(1-y(29)))*y(5)*(-params(16))/y(5)));
  g1(10,9)=1;
  g1(10,10)=(-(1-T47));
  g1(10,30)=y(10)*y(5)*(-params(16))/y(5);
  g1(11,10)=1;
  g1(11,11)=(-(params(20)*1/params(12)*getPowerDeriv(y(11)/params(12),params(3),1)));
  g1(12,4)=(-(y(28)/mu__));
  g1(12,12)=1;
  g1(12,28)=(-(y(4)/mu__));
  g1(13,2)=(-(T66*params(21)));
  g1(13,3)=(-(params(21)*y(2)*(-y(9))/(y(3)*y(3))));
  g1(13,9)=(-(params(21)*((-1)+y(2)*1/y(3))));
  g1(13,12)=(-params(21));
  g1(13,13)=1-(1-params(21))*params(2);
  g1(14,8)=(-(1-y(30)));
  g1(14,14)=1-(1-y(30));
  g1(14,30)=y(14)+y(8);
  g1(15,2)=(-(y(15)*(-(1-y(30)))));
  g1(15,15)=1-(1-y(2)*(1-y(30))-y(17));
  g1(15,17)=(-(1-y(15)));
  g1(15,30)=(-(y(15)*y(2)));
  g1(16,3)=(-((1-y(30))*(-y(16))));
  g1(16,11)=(-1);
  g1(16,15)=(-((1-y(30))*(-y(29))));
  g1(16,16)=1-(1-y(30))*(1-y(3));
  g1(16,29)=(-((1-y(30))*(1-y(15))));
  g1(16,30)=y(16)*(1-y(3))+y(29)*(1-y(15));
  g1(17,17)=1;
  g1(17,29)=(-(1-y(30)));
  g1(17,30)=(-(1-y(29)));
  g1(18,15)=1;
  g1(18,18)=1;
  g1(19,18)=1;
  g1(19,19)=(-1);
  g1(19,20)=(-1);
  g1(20,4)=(-(y(28)*y(20)));
  g1(20,20)=(-(y(4)*y(28)));
  g1(20,23)=1;
  g1(20,28)=(-(y(4)*y(20)));
  g1(21,3)=(-(y(16)*params(17)));
  g1(21,10)=(-(y(11)/(1+params(3))));
  g1(21,11)=(-(y(10)/(1+params(3))));
  g1(21,16)=(-(y(3)*params(17)));
  g1(21,21)=1;
  g1(22,21)=(-1);
  g1(22,22)=(-1);
  g1(22,23)=1;
  g1(23,6)=(-y(8));
  g1(23,8)=(-y(6));
  g1(23,23)=(-1);
  g1(23,24)=1;
  g1(24,7)=(-y(14));
  g1(24,12)=(-y(18));
  g1(24,14)=(-y(7));
  g1(24,18)=(-y(12));
  g1(24,24)=1;
  g1(25,21)=1;
  g1(25,24)=(-1);
  g1(25,25)=1;
  g1(26,13)=(-(y(18)/y(25)));
  g1(26,18)=(-(y(13)/y(25)));
  g1(26,25)=(-((-(y(13)*y(18)))/(y(25)*y(25))));
  g1(26,26)=1;
  g1(27,4)=(-((-(y(18)*y(25)))/(y(4)*y(18)*y(4)*y(18))));
  g1(27,18)=(-((-(y(4)*y(25)))/(y(4)*y(18)*y(4)*y(18))));
  g1(27,25)=(-(1/(y(4)*y(18))));
  g1(27,27)=1;
  g1(28,28)=1/y(28)-(params(5)*1/y(28)+(1-params(5))*1/(y(28)));
  g1(29,29)=1/y(29)-(params(6)*1/y(29)+(1-params(6))*1/(y(29)));
  g1(30,30)=1/y(30)-(params(7)*1/y(30)+(1-params(7))*1/(y(30)));
  g1(31,15)=(-(0.6666666666666666/(0.3333333333333333*(y(15)+y(15)+y(45)))));
  g1(31,31)=1;
  g1(31,45)=(-(0.3333333333333333/(0.3333333333333333*(y(15)+y(15)+y(45)))));
  g1(32,16)=(-(0.6666666666666666/(0.3333333333333333*(y(16)+y(16)+y(46)))));
  g1(32,32)=1;
  g1(32,46)=(-(0.3333333333333333/(0.3333333333333333*(y(16)+y(16)+y(46)))));
  g1(33,17)=(-(0.6666666666666666/(0.3333333333333333*(y(17)+y(17)+y(47)))));
  g1(33,33)=1;
  g1(33,47)=(-(0.3333333333333333/(0.3333333333333333*(y(17)+y(17)+y(47)))));
  g1(34,2)=(-(0.3333333333333333*(1-y(30)+1-y(30))/T228));
  g1(34,30)=(-(0.3333333333333333*((-y(2))-y(2))/T228));
  g1(34,34)=1;
  g1(34,48)=(-(0.3333333333333333*(-y(49))/T228));
  g1(34,49)=(-(0.3333333333333333*(1-y(48))/T228));
  g1(35,30)=(-(2/(y(48)+y(30)+y(30))));
  g1(35,35)=1;
  g1(35,48)=(-(1/(y(48)+y(30)+y(30))));
  g1(36,8)=(-(0.6666666666666666/(0.3333333333333333*(y(8)+y(8)+y(50)))));
  g1(36,36)=1;
  g1(36,50)=(-(0.3333333333333333/(0.3333333333333333*(y(8)+y(8)+y(50)))));
  g1(37,27)=(-(0.6666666666666666/(0.3333333333333333*(y(27)+y(27)+y(51)))));
  g1(37,37)=1;
  g1(37,51)=(-(0.3333333333333333/(0.3333333333333333*(y(27)+y(27)+y(51)))));
  g1(38,38)=1;
  g1(38,53)=(-1);
  g1(39,39)=1;
  g1(39,55)=(-1);
  g1(40,40)=1;
  g1(40,57)=(-1);
  g1(41,41)=1;
  g1(41,59)=(-1);
  g1(42,42)=1;
  g1(42,61)=(-1);
  g1(43,43)=1;
  g1(43,63)=(-1);
  g1(44,44)=1;
  g1(44,65)=(-1);
  g1(45,15)=(-1);
  g1(45,45)=1;
  g1(46,16)=(-1);
  g1(46,46)=1;
  g1(47,17)=(-1);
  g1(47,47)=1;
  g1(48,30)=(-1);
  g1(48,48)=1;
  g1(49,2)=(-1);
  g1(49,49)=1;
  g1(50,8)=(-1);
  g1(50,50)=1;
  g1(51,27)=(-1);
  g1(51,51)=1;
  g1(52,31)=(-1);
  g1(52,52)=1;
  g1(53,52)=(-1);
  g1(53,53)=1;
  g1(54,32)=(-1);
  g1(54,54)=1;
  g1(55,54)=(-1);
  g1(55,55)=1;
  g1(56,33)=(-1);
  g1(56,56)=1;
  g1(57,56)=(-1);
  g1(57,57)=1;
  g1(58,34)=(-1);
  g1(58,58)=1;
  g1(59,58)=(-1);
  g1(59,59)=1;
  g1(60,35)=(-1);
  g1(60,60)=1;
  g1(61,60)=(-1);
  g1(61,61)=1;
  g1(62,36)=(-1);
  g1(62,62)=1;
  g1(63,62)=(-1);
  g1(63,63)=1;
  g1(64,37)=(-1);
  g1(64,64)=1;
  g1(65,64)=(-1);
  g1(65,65)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],65,4225);
end
end
