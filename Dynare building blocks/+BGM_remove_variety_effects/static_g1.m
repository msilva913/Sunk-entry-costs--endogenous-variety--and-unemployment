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
    T = BGM_remove_variety_effects.static_g1_tt(T, y, x, params);
end
g1 = zeros(19, 19);
g1(1,1)=T(2)-(y(5)+y(4))*params(1)*(1-params(3))*T(2)/y(5);
g1(1,4)=(-(T(1)*params(1)*(1-params(3))/y(5)));
g1(1,5)=(-((T(1)*params(1)*(1-params(3))*y(5)-T(1)*params(1)*(1-params(3))*(y(5)+y(4)))/(y(5)*y(5))));
g1(2,12)=1;
g1(3,1)=1;
g1(3,8)=(-(1.0*y(18)));
g1(3,18)=(-(1.0*y(8)));
g1(4,7)=(-1);
g1(4,8)=1;
g1(4,9)=1;
g1(5,6)=1;
g1(5,18)=(-(1.0/(params(6)/(params(6)-1))));
g1(6,1)=(-(y(6)*T(2)/params(4)*T(3)));
g1(6,6)=(-(T(3)*T(1)/params(4)));
g1(6,7)=1;
g1(7,5)=1;
g1(7,6)=(-(params(7)/y(18)));
g1(7,18)=(-((-(y(6)*params(7)))/(y(18)*y(18))));
g1(8,1)=(-(1/(params(6)*y(3))));
g1(8,3)=(-((-(params(6)*y(1)))/(params(6)*y(3)*params(6)*y(3))));
g1(8,4)=1;
g1(9,1)=1;
g1(9,2)=y(5);
g1(9,3)=(-y(4));
g1(9,4)=(-y(3));
g1(9,5)=y(2);
g1(9,6)=(-y(7));
g1(9,7)=(-y(6));
g1(10,1)=(-1);
g1(10,2)=(-y(5));
g1(10,5)=(-y(2));
g1(10,10)=1;
g1(11,10)=(-1);
g1(11,11)=1;
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

end
