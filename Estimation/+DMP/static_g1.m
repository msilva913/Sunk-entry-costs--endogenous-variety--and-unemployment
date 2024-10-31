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
    T = DMP.static_g1_tt(T, y, x, params);
end
g1 = zeros(10, 10);
g1(1,1)=params(2);
g1(1,6)=(-T(2))/(y(6)*y(6))-params(2)*(-(T(2)*(1-params(4))))/(y(6)*y(6));
g1(1,7)=(-(params(2)*exp(y(7))));
g1(2,1)=1;
g1(2,4)=(-(params(7)*T(2)));
g1(2,7)=(-(params(7)*exp(y(7))));
g1(3,3)=1-(1-y(5)-params(4));
g1(3,5)=y(3);
g1(4,2)=(-(1/y(3)));
g1(4,3)=(-((-y(2))/(y(3)*y(3))));
g1(4,4)=1;
g1(5,4)=(-(T(1)*getPowerDeriv(y(4),1-params(9),1)));
g1(5,5)=1;
g1(6,4)=(-(T(1)*getPowerDeriv(y(4),(-params(9)),1)));
g1(6,6)=1;
g1(7,3)=(-(100*((y(3))-y(3))/((y(3))*(y(3)))/(y(3)/(y(3)))));
g1(7,8)=1;
g1(8,2)=(-(100*((y(2))-y(2))/((y(2))*(y(2)))/(y(2)/(y(2)))));
g1(8,9)=1;
g1(9,4)=(-(100*((y(4))-y(4))/((y(4))*(y(4)))/(y(4)/(y(4)))));
g1(9,10)=1;
g1(10,7)=1-params(10);
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end
