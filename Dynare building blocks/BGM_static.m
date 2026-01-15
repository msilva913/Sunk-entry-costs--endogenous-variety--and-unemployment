function [residual, g1, g2] = BGM_static(y, x, params)
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

residual = zeros( 19, 1);

%
% Model equations
%

mu__ = params(6)/(params(6)-1);
T14 = y(1)^(-params(8));
T57 = 1/y(6);
T103 = getPowerDeriv(y(1),(-params(8)),1);
T110 = getPowerDeriv(T14*y(6)/params(4),params(5),1);
lhs =T14;
rhs =T14*params(1)*(1-params(3))*(y(5)+y(4))/y(5);
residual(1)= lhs-rhs;
lhs =y(12);
rhs =y(3)^(1/(params(6)-1));
residual(2)= lhs-rhs;
lhs =y(1);
rhs =y(12)*y(18)*y(8);
residual(3)= lhs-rhs;
lhs =y(9);
rhs =y(7)-y(8);
residual(4)= lhs-rhs;
lhs =y(6);
rhs =y(12)*y(18)/mu__;
residual(5)= lhs-rhs;
lhs =y(7);
rhs =(T14*y(6)/params(4))^params(5);
residual(6)= lhs-rhs;
lhs =y(5);
rhs =y(6)*params(7)/y(18);
residual(7)= lhs-rhs;
lhs =y(4);
rhs =y(1)/(params(6)*y(3));
residual(8)= lhs-rhs;
lhs =y(7);
rhs =T57*(y(1)/mu__+y(5)*y(2));
residual(9)= lhs-rhs;
lhs =y(10);
rhs =y(1)+y(5)*y(2);
residual(10)= lhs-rhs;
lhs =y(11);
rhs =y(10)/y(12);
residual(11)= lhs-rhs;
lhs =y(3);
rhs =(1-params(3))*(y(3)+y(2));
residual(12)= lhs-rhs;
lhs =log(y(18));
rhs =log(y(18))*params(9)+x(1);
residual(13)= lhs-rhs;
lhs =y(13);
rhs =100*log(y(10));
residual(14)= lhs-rhs;
lhs =y(14);
rhs =100*log(y(11));
residual(15)= lhs-rhs;
lhs =y(15);
rhs =100*log(y(7));
residual(16)= lhs-rhs;
lhs =y(16);
rhs =100*log(y(8));
residual(17)= lhs-rhs;
lhs =y(17);
rhs =100*log(y(9));
residual(18)= lhs-rhs;
lhs =y(19);
rhs =log(y(18))*100;
residual(19)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(19, 19);

  %
  % Jacobian matrix
  %

  g1(1,1)=T103-(y(5)+y(4))*params(1)*(1-params(3))*T103/y(5);
  g1(1,4)=(-(T14*params(1)*(1-params(3))/y(5)));
  g1(1,5)=(-((T14*params(1)*(1-params(3))*y(5)-T14*params(1)*(1-params(3))*(y(5)+y(4)))/(y(5)*y(5))));
  g1(2,3)=(-(getPowerDeriv(y(3),1/(params(6)-1),1)));
  g1(2,12)=1;
  g1(3,1)=1;
  g1(3,8)=(-(y(12)*y(18)));
  g1(3,12)=(-(y(18)*y(8)));
  g1(3,18)=(-(y(12)*y(8)));
  g1(4,7)=(-1);
  g1(4,8)=1;
  g1(4,9)=1;
  g1(5,6)=1;
  g1(5,12)=(-(y(18)/mu__));
  g1(5,18)=(-(y(12)/mu__));
  g1(6,1)=(-(y(6)*T103/params(4)*T110));
  g1(6,6)=(-(T110*T14/params(4)));
  g1(6,7)=1;
  g1(7,5)=1;
  g1(7,6)=(-(params(7)/y(18)));
  g1(7,18)=(-((-(y(6)*params(7)))/(y(18)*y(18))));
  g1(8,1)=(-(1/(params(6)*y(3))));
  g1(8,3)=(-((-(params(6)*y(1)))/(params(6)*y(3)*params(6)*y(3))));
  g1(8,4)=1;
  g1(9,1)=(-(T57*1/mu__));
  g1(9,2)=(-(y(5)*T57));
  g1(9,5)=(-(T57*y(2)));
  g1(9,6)=(-((y(1)/mu__+y(5)*y(2))*(-1)/(y(6)*y(6))));
  g1(9,7)=1;
  g1(10,1)=(-1);
  g1(10,2)=(-y(5));
  g1(10,5)=(-y(2));
  g1(10,10)=1;
  g1(11,10)=(-(1/y(12)));
  g1(11,11)=1;
  g1(11,12)=(-((-y(10))/(y(12)*y(12))));
  g1(12,2)=(-(1-params(3)));
  g1(12,3)=1-(1-params(3));
  g1(13,18)=1/y(18)-params(9)*1/y(18);
  g1(14,10)=(-(100*1/y(10)));
  g1(14,13)=1;
  g1(15,11)=(-(100*1/y(11)));
  g1(15,14)=1;
  g1(16,7)=(-(100*1/y(7)));
  g1(16,15)=1;
  g1(17,8)=(-(100*1/y(8)));
  g1(17,16)=1;
  g1(18,9)=(-(100*1/y(9)));
  g1(18,17)=1;
  g1(19,18)=(-(100*1/y(18)));
  g1(19,19)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],19,361);
end
end
