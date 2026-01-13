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
    T = BGM_translog.static_g1_tt(T, y, x, params);
end
g1 = zeros(21, 21);
g1(1,3)=(-((-params(9))/(params(9)*y(3)*params(9)*y(3))));
g1(1,20)=1;
g1(2,20)=(-((y(20)-1-y(20))/((y(20)-1)*(y(20)-1))));
g1(2,21)=1;
g1(3,3)=(-(T(1)*1/(sqrt(y(3))+sqrt(y(3)))+sqrt(y(3))*T(1)*(-(2*params(9)))/(y(3)*2*params(9)*y(3)*2*params(9))));
g1(3,12)=1;
g1(4,1)=T(3)-(y(5)+y(4))*params(1)*(1-params(3))*T(3)/y(5);
g1(4,4)=(-(T(2)*params(1)*(1-params(3))/y(5)));
g1(4,5)=(-((T(2)*params(1)*(1-params(3))*y(5)-T(2)*params(1)*(1-params(3))*(y(5)+y(4)))/(y(5)*y(5))));
g1(5,1)=1;
g1(5,8)=(-(y(12)*y(18)));
g1(5,12)=(-(y(18)*y(8)));
g1(5,18)=(-(y(12)*y(8)));
g1(6,6)=1;
g1(6,12)=(-(y(18)/y(20)));
g1(6,18)=(-(y(12)/y(20)));
g1(6,20)=(-((-(y(12)*y(18)))/(y(20)*y(20))));
g1(7,1)=(-(y(6)*T(3)/params(4)*T(4)));
g1(7,6)=(-(T(4)*T(2)/params(4)));
g1(7,7)=1;
g1(8,7)=(-1);
g1(8,8)=1;
g1(8,9)=1;
g1(9,5)=1;
g1(9,6)=(-(params(6)/y(18)));
g1(9,18)=(-((-(y(6)*params(6)))/(y(18)*y(18))));
g1(10,1)=(-((y(20)-1)/y(20)*1/y(3)));
g1(10,3)=(-((y(20)-1)/y(20)*(-y(1))/(y(3)*y(3))));
g1(10,4)=1;
g1(10,20)=(-(y(1)/y(3)*(y(20)-(y(20)-1))/(y(20)*y(20))));
g1(11,1)=1;
g1(11,2)=y(5);
g1(11,3)=(-y(4));
g1(11,4)=(-y(3));
g1(11,5)=y(2);
g1(11,6)=(-y(7));
g1(11,7)=(-y(6));
g1(12,1)=(-1);
g1(12,2)=(-y(5));
g1(12,5)=(-y(2));
g1(12,10)=1;
g1(13,10)=(-(1/y(12)));
g1(13,11)=1;
g1(13,12)=(-((-y(10))/(y(12)*y(12))));
g1(14,2)=(-(1-params(3)));
g1(14,3)=1-(1-params(3));
g1(15,18)=1/y(18)-params(8)*1/y(18);
g1(16,10)=(-(100*1/y(10)));
g1(16,13)=1;
g1(17,11)=(-(100*1/y(11)));
g1(17,14)=1;
g1(18,7)=(-(100*1/y(7)));
g1(18,15)=1;
g1(19,8)=(-(100*1/y(8)));
g1(19,16)=1;
g1(20,9)=(-(100*1/y(9)));
g1(20,17)=1;
g1(21,18)=(-(100*1/y(18)));
g1(21,19)=1;

end
