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
    T = baseline_alt.static_resid_tt(T, y, x, params);
end
residual = zeros(37, 1);
lhs = y(2)/y(9);
rhs = params(2)*(1-params(3))*(y(13)-y(3)-y(2)+y(2)*(1-(params(6)-params(3))/(1-params(3)))/y(9));
residual(1) = lhs - rhs;
lhs = y(13);
rhs = T(16)/T(2);
residual(2) = lhs - rhs;
lhs = y(3);
rhs = params(10)*(y(13)-y(2)+y(2)*y(7)/(1-params(3)))+params(11)*params(1)*(1-params(10));
residual(3) = lhs - rhs;
lhs = y(1);
rhs = T(17)^params(9);
residual(4) = lhs - rhs;
lhs = y(2);
rhs = y(1)-params(2)*(1-params(3))*y(1);
residual(5) = lhs - rhs;
lhs = y(7);
rhs = y(5)/y(6);
residual(6) = lhs - rhs;
lhs = y(8);
rhs = T(4)*y(7)^(1-params(12));
residual(7) = lhs - rhs;
lhs = y(9);
rhs = T(4)*y(7)^(-params(12));
residual(8) = lhs - rhs;
lhs = y(21);
rhs = 1-y(6);
residual(9) = lhs - rhs;
lhs = y(21);
rhs = y(23)+y(22);
residual(10) = lhs - rhs;
lhs = y(12);
rhs = y(10)^(1/(params(4)-1));
residual(11) = lhs - rhs;
lhs = y(19);
rhs = y(17)^(-params(5));
residual(12) = lhs - rhs;
lhs = y(14);
rhs = params(2)*(1-params(3))*(y(14)+y(15));
residual(13) = lhs - rhs;
lhs = y(16);
rhs = T(16)*y(23);
residual(14) = lhs - rhs;
lhs = y(16);
rhs = y(17)+y(18);
residual(15) = lhs - rhs;
lhs = y(18);
rhs = T(18)*T(19);
residual(16) = lhs - rhs;
lhs = y(11);
rhs = exp(y(25))*T(12)*y(22)/T(13);
residual(17) = lhs - rhs;
lhs = y(14);
rhs = T(13)*y(12)/T(2);
residual(18) = lhs - rhs;
lhs = y(20);
rhs = y(16)+y(14)*y(11);
residual(19) = lhs - rhs;
lhs = y(20);
rhs = y(13)*y(21)+y(10)*y(15);
residual(20) = lhs - rhs;
lhs = y(5);
rhs = y(4)+(1-params(3))*(y(5)*(1-y(9))+(params(6)-params(3))/(1-params(3))*(1-y(6)));
residual(21) = lhs - rhs;
lhs = y(6);
rhs = y(6)*(1-(1-params(3))*y(8))+params(6)*(1-y(6));
residual(22) = lhs - rhs;
lhs = y(10);
rhs = (1-params(3))*(y(10)+y(11));
residual(23) = lhs - rhs;
lhs = y(24);
rhs = y(3)*y(21)/y(20);
residual(24) = lhs - rhs;
residual(25) = y(26);
lhs = y(27);
rhs = y(26)+y(26)+y(26);
residual(26) = lhs - rhs;
lhs = y(28);
rhs = y(6)-(y(6));
residual(27) = lhs - rhs;
lhs = y(29);
rhs = y(5)-(y(5));
residual(28) = lhs - rhs;
lhs = y(30);
rhs = 100*log(y(6));
residual(29) = lhs - rhs;
lhs = y(31);
rhs = 100*log(y(5));
residual(30) = lhs - rhs;
lhs = y(32);
rhs = 100*log(y(7));
residual(31) = lhs - rhs;
lhs = y(33);
rhs = 100*log(y(4));
residual(32) = lhs - rhs;
lhs = y(34);
rhs = 100*log(y(10));
residual(33) = lhs - rhs;
lhs = y(35);
rhs = 100*log(y(11));
residual(34) = lhs - rhs;
lhs = y(36);
rhs = 100*log(y(17));
residual(35) = lhs - rhs;
lhs = y(25);
rhs = y(25)*params(14)-x(1);
residual(36) = lhs - rhs;
lhs = y(37);
rhs = y(26);
residual(37) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
