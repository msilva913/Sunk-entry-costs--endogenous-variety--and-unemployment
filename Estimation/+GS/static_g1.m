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
    T = GS.static_g1_tt(T, y, x, params);
end
g1 = zeros(23, 23);
g1(1,1)=1/y(8)-T(6)*(T(7)/y(8)-1);
g1(1,2)=T(6);
g1(1,8)=(-y(1))/(y(8)*y(8))-T(6)*(-(y(1)*T(7)))/(y(8)*y(8));
g1(1,10)=(-(T(6)*exp(y(10))));
g1(1,11)=(-(T(6)*(exp(y(10))-y(2)-y(1)+y(1)*T(7)/y(8))));
g1(1,13)=(-((exp(y(10))-y(2)-y(1)+y(1)*T(7)/y(8))*params(2)*exp(y(11))*(-(params(3)*exp(y(13))))));
g1(2,1)=1;
g1(2,3)=(-((T(11)-T(6)*T(11))/T(9)));
g1(2,11)=(-((-(T(6)*T(8)))/T(9)));
g1(2,13)=(-((-(T(8)*params(2)*exp(y(11))*(-(params(3)*exp(y(13))))))/T(9)));
g1(3,1)=(-(params(8)*exp(y(12))*(y(6)/T(5)-1)));
g1(3,2)=1;
g1(3,6)=(-(params(8)*exp(y(12))*y(1)*1/T(5)));
g1(3,10)=(-(exp(y(10))*params(8)*exp(y(12))));
g1(3,12)=(-(T(10)+params(9)*(-(params(8)*exp(y(12))))));
g1(3,13)=(-(params(8)*exp(y(12))*y(1)*(-(y(6)*(-(params(3)*exp(y(13))))))/(T(5)*T(5))));
g1(4,3)=(-1);
g1(4,4)=1-T(5)*(1-y(8));
g1(4,5)=(-(T(5)*(-T(4))));
g1(4,8)=(-(T(5)*(-y(4))));
g1(4,13)=(-((y(4)*(1-y(8))+T(4)*(1-y(5)))*(-(params(3)*exp(y(13))))));
g1(5,5)=1-(1-T(5)*y(7)-y(9));
g1(5,7)=(-(y(5)*(-T(5))));
g1(5,9)=(-(1-y(5)));
g1(5,13)=(-(y(5)*(-(y(7)*(-(params(3)*exp(y(13))))))));
g1(6,4)=(-(1/y(5)));
g1(6,5)=(-((-y(4))/(y(5)*y(5))));
g1(6,6)=1;
g1(7,6)=(-(T(3)*getPowerDeriv(y(6),1-params(10),1)));
g1(7,7)=1;
g1(8,6)=(-(T(3)*getPowerDeriv(y(6),(-params(10)),1)));
g1(8,8)=1;
g1(9,9)=1;
g1(9,13)=T(7)*(-(params(3)*exp(y(13))));
g1(10,14)=1;
g1(11,14)=(-3);
g1(11,15)=1;
g1(12,16)=1;
g1(13,17)=1;
g1(14,9)=(-1);
g1(14,18)=1;
g1(15,5)=(-(100*((y(5))-y(5))/((y(5))*(y(5)))/(y(5)/(y(5)))));
g1(15,19)=1;
g1(16,4)=(-(100*((y(4))-y(4))/((y(4))*(y(4)))/(y(4)/(y(4)))));
g1(16,20)=1;
g1(17,6)=(-(100*((y(6))-y(6))/((y(6))*(y(6)))/(y(6)/(y(6)))));
g1(17,21)=1;
g1(18,3)=(-(100*((y(3))-y(3))/((y(3))*(y(3)))/(y(3)/(y(3)))));
g1(18,22)=1;
g1(19,10)=1-params(11);
g1(20,11)=1-params(12);
g1(21,12)=1-params(13);
g1(22,13)=1-params(14);
g1(23,14)=(-1);
g1(23,23)=1;
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end
