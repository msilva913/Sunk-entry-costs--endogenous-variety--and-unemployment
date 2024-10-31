function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
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
%   residual
%

if T_flag
    T = baseline_u_growth_est.static_resid_tt(T, y, x, params);
end
residual = zeros(49, 1);
lhs = y(2)/y(9);
rhs = T(14);
residual(1) = lhs - rhs;
lhs = y(13);
rhs = y(12)*exp(y(28))/T(4);
residual(2) = lhs - rhs;
lhs = y(3);
rhs = T(15)*T(16)+T(7)*params(10)*(1-T(16));
residual(3) = lhs - rhs;
lhs = y(1);
rhs = T(17)^params(9);
residual(4) = lhs - rhs;
lhs = y(2);
rhs = y(1)-(1-params(3))*params(2)*exp(y(29))*y(1);
residual(5) = lhs - rhs;
lhs = y(7);
rhs = y(5)/y(6);
residual(6) = lhs - rhs;
lhs = y(8);
rhs = T(19);
residual(7) = lhs - rhs;
lhs = y(9);
rhs = T(20);
residual(8) = lhs - rhs;
lhs = y(22);
rhs = 1-y(6);
residual(9) = lhs - rhs;
lhs = y(22);
rhs = y(24)+y(23);
residual(10) = lhs - rhs;
lhs = y(12);
rhs = y(10)^T(5);
residual(11) = lhs - rhs;
lhs = y(19);
rhs = y(17)^(-params(5));
residual(12) = lhs - rhs;
lhs = y(14);
rhs = (1-params(3))*params(2)*exp(y(29))*(y(14)+y(15));
residual(13) = lhs - rhs;
lhs = y(16);
rhs = y(12)*exp(y(28))*y(24);
residual(14) = lhs - rhs;
lhs = y(16);
rhs = y(17)+y(18);
residual(15) = lhs - rhs;
lhs = y(18);
rhs = T(21)*T(22);
residual(16) = lhs - rhs;
lhs = y(11);
rhs = exp(y(28))*y(23)/T(13);
residual(17) = lhs - rhs;
lhs = y(14);
rhs = y(12)*T(13)/T(4);
residual(18) = lhs - rhs;
lhs = y(16)+y(14)*y(11);
rhs = y(13)*y(22)+y(10)*y(15);
residual(19) = lhs - rhs;
lhs = y(20);
rhs = y(16)+y(14)*y(11);
residual(20) = lhs - rhs;
lhs = y(21);
rhs = y(20)/y(12);
residual(21) = lhs - rhs;
lhs = y(5);
rhs = y(4)+(1-params(3))*(y(5)*(1-y(9))+(params(6)-params(3))/(1-params(3))*(1-y(6)));
residual(22) = lhs - rhs;
lhs = y(6);
rhs = y(6)*(1-(1-params(3))*y(8))+params(6)*(1-y(6));
residual(23) = lhs - rhs;
lhs = y(10);
rhs = (1-params(3))*(y(10)+y(11));
residual(24) = lhs - rhs;
lhs = y(33);
rhs = y(21)/y(22);
residual(25) = lhs - rhs;
lhs = y(25);
rhs = y(3)*y(22)/y(20);
residual(26) = lhs - rhs;
lhs = y(26);
rhs = y(17)/y(20);
residual(27) = lhs - rhs;
lhs = y(27);
rhs = (y(26)+y(26)+y(26))/3;
residual(28) = lhs - rhs;
residual(29) = y(34);
residual(30) = y(36);
lhs = y(35);
rhs = y(34)+y(34)+y(34);
residual(31) = lhs - rhs;
lhs = y(37);
rhs = y(36)+y(36)+y(36);
residual(32) = lhs - rhs;
residual(33) = y(38);
residual(34) = y(39);
lhs = y(40);
rhs = 100*log(y(6));
residual(35) = lhs - rhs;
lhs = y(41);
rhs = 100*log(y(5));
residual(36) = lhs - rhs;
lhs = y(42);
rhs = 100*log(y(7));
residual(37) = lhs - rhs;
lhs = y(43);
rhs = 100*log(y(4));
residual(38) = lhs - rhs;
lhs = y(44);
rhs = 100*log(y(10));
residual(39) = lhs - rhs;
lhs = y(45);
rhs = 100*log(y(11));
residual(40) = lhs - rhs;
lhs = y(46);
rhs = 100*log(y(17));
residual(41) = lhs - rhs;
lhs = y(28);
rhs = y(28)*params(14)-x(1);
residual(42) = lhs - rhs;
lhs = y(29);
rhs = y(29)*params(15)+x(2);
residual(43) = lhs - rhs;
lhs = y(30);
rhs = y(30)*params(16)+x(3);
residual(44) = lhs - rhs;
lhs = y(31);
rhs = y(31)*params(17)+x(4);
residual(45) = lhs - rhs;
lhs = y(32);
rhs = y(32)*params(18)+x(5);
residual(46) = lhs - rhs;
lhs = y(47);
rhs = y(26);
residual(47) = lhs - rhs;
lhs = y(48);
rhs = y(34);
residual(48) = lhs - rhs;
lhs = y(49);
rhs = y(36);
residual(49) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
