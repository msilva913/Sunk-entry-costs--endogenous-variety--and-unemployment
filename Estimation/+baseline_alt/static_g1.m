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
    T = baseline_alt.static_g1_tt(T, y, x, params);
end
g1 = zeros(37, 37);
g1(1,2)=1/y(9)-params(2)*(1-params(3))*((1-(params(6)-params(3))/(1-params(3)))/y(9)-1);
g1(1,3)=params(2)*(1-params(3));
g1(1,9)=(-y(2))/(y(9)*y(9))-params(2)*(1-params(3))*(-(y(2)*(1-(params(6)-params(3))/(1-params(3)))))/(y(9)*y(9));
g1(1,13)=(-(params(2)*(1-params(3))));
g1(2,7)=(-(exp(y(25))*y(12)*T(11)*(T(23)+T(24))/T(2)));
g1(2,9)=(-(T(32)/T(2)));
g1(2,12)=(-(T(12)*exp(y(25))/T(2)));
g1(2,13)=1;
g1(2,25)=(-(T(16)/T(2)));
g1(3,2)=(-(params(10)*(y(7)/(1-params(3))-1)));
g1(3,3)=1;
g1(3,7)=(-(params(10)*y(2)*1/(1-params(3))));
g1(3,13)=(-params(10));
g1(4,1)=1;
g1(4,4)=(-(1/T(15)*T(20)));
g1(4,7)=(-(T(20)*(-(y(4)*T(26)))/(T(15)*T(15))));
g1(4,9)=(-(T(20)*(-(y(4)*T(33)))/(T(15)*T(15))));
g1(5,1)=(-(1-params(2)*(1-params(3))));
g1(5,2)=1;
g1(6,5)=(-(1/y(6)));
g1(6,6)=(-((-y(5))/(y(6)*y(6))));
g1(6,7)=1;
g1(7,7)=(-(T(4)*getPowerDeriv(y(7),1-params(12),1)));
g1(7,8)=1;
g1(8,7)=(-(T(4)*getPowerDeriv(y(7),(-params(12)),1)));
g1(8,9)=1;
g1(9,6)=1;
g1(9,21)=1;
g1(10,21)=1;
g1(10,22)=(-1);
g1(10,23)=(-1);
g1(11,10)=(-(getPowerDeriv(y(10),1/(params(4)-1),1)));
g1(11,12)=1;
g1(12,17)=(-(getPowerDeriv(y(17),(-params(5)),1)));
g1(12,19)=1;
g1(13,14)=1-params(2)*(1-params(3));
g1(13,15)=(-(params(2)*(1-params(3))));
g1(14,7)=(-(y(23)*exp(y(25))*y(12)*T(11)*(T(23)+T(24))));
g1(14,9)=(-(y(23)*T(32)));
g1(14,12)=(-(y(23)*T(12)*exp(y(25))));
g1(14,16)=1;
g1(14,23)=(-T(16));
g1(14,25)=(-(T(16)*y(23)));
g1(15,16)=1;
g1(15,17)=(-1);
g1(15,18)=(-1);
g1(16,4)=(-(T(18)*1/T(15)*T(21)));
g1(16,7)=(-(T(19)*T(26)/(1+params(9))+T(18)*T(21)*(-(y(4)*T(26)))/(T(15)*T(15))));
g1(16,9)=(-(T(19)*T(33)/(1+params(9))+T(18)*T(21)*(-(y(4)*T(33)))/(T(15)*T(15))));
g1(16,18)=1;
g1(17,7)=(-((T(13)*exp(y(25))*y(22)*T(11)*(T(23)+T(24))-exp(y(25))*T(12)*y(22)*T(27))/(T(13)*T(13))));
g1(17,9)=(-((T(13)*exp(y(25))*y(22)*T(31)-exp(y(25))*T(12)*y(22)*T(34))/(T(13)*T(13))));
g1(17,11)=1;
g1(17,22)=(-(T(12)*exp(y(25))/T(13)));
g1(17,25)=(-(exp(y(25))*T(12)*y(22)/T(13)));
g1(18,7)=(-(y(12)*T(27)/T(2)));
g1(18,9)=(-(y(12)*T(34)/T(2)));
g1(18,12)=(-(T(13)/T(2)));
g1(18,14)=1;
g1(19,11)=(-y(14));
g1(19,14)=(-y(11));
g1(19,16)=(-1);
g1(19,20)=1;
g1(20,10)=(-y(15));
g1(20,13)=(-y(21));
g1(20,15)=(-y(10));
g1(20,20)=1;
g1(20,21)=(-y(13));
g1(21,4)=(-1);
g1(21,5)=1-(1-params(3))*(1-y(9));
g1(21,6)=(-((1-params(3))*(-((params(6)-params(3))/(1-params(3))))));
g1(21,9)=(-((1-params(3))*(-y(5))));
g1(22,6)=1-(1-(1-params(3))*y(8)-params(6));
g1(22,8)=(-(y(6)*(-(1-params(3)))));
g1(23,10)=1-(1-params(3));
g1(23,11)=(-(1-params(3)));
g1(24,3)=(-(y(21)/y(20)));
g1(24,20)=(-((-(y(3)*y(21)))/(y(20)*y(20))));
g1(24,21)=(-(y(3)/y(20)));
g1(24,24)=1;
g1(25,26)=1;
g1(26,26)=(-3);
g1(26,27)=1;
g1(27,28)=1;
g1(28,29)=1;
g1(29,6)=(-(100*1/y(6)));
g1(29,30)=1;
g1(30,5)=(-(100*1/y(5)));
g1(30,31)=1;
g1(31,7)=(-(100*1/y(7)));
g1(31,32)=1;
g1(32,4)=(-(100*1/y(4)));
g1(32,33)=1;
g1(33,10)=(-(100*1/y(10)));
g1(33,34)=1;
g1(34,11)=(-(100*1/y(11)));
g1(34,35)=1;
g1(35,17)=(-(100*1/y(17)));
g1(35,36)=1;
g1(36,25)=1-params(14);
g1(37,26)=(-1);
g1(37,37)=1;
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end
