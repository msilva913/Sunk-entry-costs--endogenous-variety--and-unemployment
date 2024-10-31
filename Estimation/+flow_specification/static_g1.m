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
    T = flow_specification.static_g1_tt(T, y, x, params);
end
g1 = zeros(32, 32);
g1(1,1)=params(2)*(1-params(3));
g1(1,6)=(-T(7))/(y(6)*y(6))-params(2)*(1-params(3))*(-(T(7)*(1-(params(6)-params(3))/(1-params(3)))))/(y(6)*y(6));
g1(1,10)=(-(params(2)*(1-params(3))));
g1(2,9)=(-(T(6)*exp(y(21))/T(3)));
g1(2,10)=1;
g1(2,21)=(-(T(10)/T(3)));
g1(3,1)=1;
g1(3,4)=(-(T(8)*T(7)/(1-params(3))));
g1(3,10)=(-T(8));
g1(4,2)=(-(1/y(3)));
g1(4,3)=(-((-y(2))/(y(3)*y(3))));
g1(4,4)=1;
g1(5,4)=(-(T(2)*getPowerDeriv(y(4),1-params(10),1)));
g1(5,5)=1;
g1(6,4)=(-(T(2)*getPowerDeriv(y(4),(-params(10)),1)));
g1(6,6)=1;
g1(7,3)=1;
g1(7,17)=1;
g1(8,17)=1;
g1(8,18)=(-1);
g1(8,19)=(-1);
g1(9,7)=(-(getPowerDeriv(y(7),T(4),1)));
g1(9,9)=1;
g1(10,14)=(-(getPowerDeriv(y(14),(-params(5)),1)));
g1(10,15)=1;
g1(11,11)=1-params(2)*(1-params(3));
g1(11,12)=(-(params(2)*(1-params(3))));
g1(12,9)=(-(y(19)*T(6)*exp(y(21))));
g1(12,13)=1;
g1(12,19)=(-T(10));
g1(12,21)=(-(T(10)*y(19)));
g1(13,2)=(-T(7));
g1(13,13)=1;
g1(13,14)=(-1);
g1(14,8)=1;
g1(14,18)=(-(T(6)*exp(y(21))/T(9)));
g1(14,21)=(-(exp(y(21))*T(6)*y(18)/T(9)));
g1(15,9)=(-(T(9)/T(3)));
g1(15,11)=1;
g1(16,8)=(-y(11));
g1(16,11)=(-y(8));
g1(16,13)=(-1);
g1(16,16)=1;
g1(17,7)=(-y(12));
g1(17,10)=(-y(17));
g1(17,12)=(-y(7));
g1(17,16)=1;
g1(17,17)=(-y(10));
g1(18,3)=1-(1-(1-params(3))*y(5)-params(6));
g1(18,5)=(-(y(3)*(-(1-params(3)))));
g1(19,7)=1-(1-params(3));
g1(19,8)=(-(1-params(3)));
g1(20,1)=(-(y(17)/y(16)));
g1(20,16)=(-((-(y(1)*y(17)))/(y(16)*y(16))));
g1(20,17)=(-(y(1)/y(16)));
g1(20,20)=1;
g1(21,22)=1;
g1(22,22)=(-3);
g1(22,23)=1;
g1(23,24)=1;
g1(24,25)=1;
g1(25,3)=(-(100*1/y(3)));
g1(25,26)=1;
g1(26,2)=(-(100*1/y(2)));
g1(26,27)=1;
g1(27,4)=(-(100*1/y(4)));
g1(27,28)=1;
g1(28,7)=(-(100*1/y(7)));
g1(28,29)=1;
g1(29,8)=(-(100*1/y(8)));
g1(29,30)=1;
g1(30,14)=(-(100*1/y(14)));
g1(30,31)=1;
g1(31,21)=1-params(13);
g1(32,22)=(-1);
g1(32,32)=1;
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end
