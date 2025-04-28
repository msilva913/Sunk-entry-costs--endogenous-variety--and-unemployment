function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g1
%

if T_flag
    T = GS_basic.static_g1_tt(T, y, x, params);
end
g1 = zeros(23, 23);
g1(1,1)=1/y(8)-(1-params(3))*params(2)*exp(y(10))*((1-T(4))/y(8)-1);
g1(1,2)=(1-params(3))*params(2)*exp(y(10));
g1(1,8)=(-y(1))/(y(8)*y(8))-(1-params(3))*params(2)*exp(y(10))*(-(y(1)*(1-T(4))))/(y(8)*y(8));
g1(1,9)=(-((1-params(3))*params(2)*exp(y(10))*exp(y(9))));
g1(1,10)=(-T(5));
g1(2,1)=1;
g1(2,3)=(-((T(9)-(1-params(3))*params(2)*exp(y(10))*T(9))/T(7)));
g1(2,10)=(-((-((1-params(3))*params(2)*exp(y(10))*T(6)))/T(7)));
g1(3,1)=(-(params(8)*exp(y(11))*(y(6)/(1-params(3))-1)));
g1(3,2)=1;
g1(3,6)=(-(params(8)*exp(y(11))*y(1)*1/(1-params(3))));
g1(3,9)=(-(exp(y(9))*params(8)*exp(y(11))));
g1(3,11)=(-(T(8)+params(9)*(-(params(8)*exp(y(11))))));
g1(4,3)=(-1);
g1(4,4)=1-(1-params(3))*(1-y(8));
g1(4,5)=(-((1-params(3))*(-T(4))));
g1(4,8)=(-((1-params(3))*(-y(4))));
g1(5,5)=1-(1-(1-params(3))*y(7)-params(4));
g1(5,7)=(-(y(5)*(-(1-params(3)))));
g1(6,4)=(-(1/y(5)));
g1(6,5)=(-((-y(4))/(y(5)*y(5))));
g1(6,6)=1;
g1(7,6)=(-(T(3)*getPowerDeriv(y(6),1-params(10),1)));
g1(7,7)=1;
g1(8,6)=(-(T(3)*getPowerDeriv(y(6),(-params(10)),1)));
g1(8,8)=1;
g1(9,12)=1;
g1(10,12)=(-2);
g1(10,13)=1;
g1(10,22)=(-1);
g1(11,14)=1;
g1(12,14)=(-2);
g1(12,15)=1;
g1(12,23)=(-1);
g1(13,16)=1;
g1(14,17)=1;
g1(15,5)=(-(100*((y(5))-y(5))/((y(5))*(y(5)))/(y(5)/(y(5)))));
g1(15,18)=1;
g1(16,4)=(-(100*((y(4))-y(4))/((y(4))*(y(4)))/(y(4)/(y(4)))));
g1(16,19)=1;
g1(17,6)=(-(100*((y(6))-y(6))/((y(6))*(y(6)))/(y(6)/(y(6)))));
g1(17,20)=1;
g1(18,3)=(-(100*((y(3))-y(3))/((y(3))*(y(3)))/(y(3)/(y(3)))));
g1(18,21)=1;
g1(19,9)=1-params(11);
g1(20,10)=1-params(12);
g1(21,11)=1-params(13);
g1(22,12)=(-1);
g1(22,22)=1;
g1(23,14)=(-1);
g1(23,23)=1;

end
