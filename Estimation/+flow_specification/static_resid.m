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
    T = flow_specification.static_resid_tt(T, y, x, params);
end
residual = zeros(32, 1);
lhs = T(7)/y(6);
rhs = params(2)*(1-params(3))*(y(10)-y(1)+T(7)*(1-(params(6)-params(3))/(1-params(3)))/y(6));
residual(1) = lhs - rhs;
lhs = y(10);
rhs = T(10)/T(3);
residual(2) = lhs - rhs;
lhs = y(1);
rhs = T(8)*(y(10)+T(7)*y(4)/(1-params(3)))+params(9)*params(1)*(1-T(8));
residual(3) = lhs - rhs;
lhs = y(4);
rhs = y(2)/y(3);
residual(4) = lhs - rhs;
lhs = y(5);
rhs = T(2)*y(4)^(1-params(10));
residual(5) = lhs - rhs;
lhs = y(6);
rhs = T(2)*y(4)^(-params(10));
residual(6) = lhs - rhs;
lhs = y(17);
rhs = 1-y(3);
residual(7) = lhs - rhs;
lhs = y(17);
rhs = y(19)+y(18);
residual(8) = lhs - rhs;
lhs = y(9);
rhs = y(7)^T(4);
residual(9) = lhs - rhs;
lhs = y(15);
rhs = y(14)^(-params(5));
residual(10) = lhs - rhs;
lhs = y(11);
rhs = params(2)*(1-params(3))*(y(11)+y(12));
residual(11) = lhs - rhs;
lhs = y(13);
rhs = T(10)*y(19);
residual(12) = lhs - rhs;
lhs = y(13);
rhs = y(14)+T(7)*y(2);
residual(13) = lhs - rhs;
lhs = y(8);
rhs = exp(y(21))*T(6)*y(18)/T(9);
residual(14) = lhs - rhs;
lhs = y(11);
rhs = T(9)*y(9)/T(3);
residual(15) = lhs - rhs;
lhs = y(16);
rhs = y(13)+y(11)*y(8);
residual(16) = lhs - rhs;
lhs = y(16);
rhs = y(10)*y(17)+y(7)*y(12);
residual(17) = lhs - rhs;
lhs = y(3);
rhs = y(3)*(1-(1-params(3))*y(5))+params(6)*(1-y(3));
residual(18) = lhs - rhs;
lhs = y(7);
rhs = (1-params(3))*(y(7)+y(8));
residual(19) = lhs - rhs;
lhs = y(20);
rhs = y(1)*y(17)/y(16);
residual(20) = lhs - rhs;
residual(21) = y(22);
lhs = y(23);
rhs = y(22)+y(22)+y(22);
residual(22) = lhs - rhs;
lhs = y(24);
rhs = y(3)-(y(3));
residual(23) = lhs - rhs;
lhs = y(25);
rhs = y(2)-(y(2));
residual(24) = lhs - rhs;
lhs = y(26);
rhs = 100*log(y(3));
residual(25) = lhs - rhs;
lhs = y(27);
rhs = 100*log(y(2));
residual(26) = lhs - rhs;
lhs = y(28);
rhs = 100*log(y(4));
residual(27) = lhs - rhs;
lhs = y(29);
rhs = 100*log(y(7));
residual(28) = lhs - rhs;
lhs = y(30);
rhs = 100*log(y(8));
residual(29) = lhs - rhs;
lhs = y(31);
rhs = 100*log(y(14));
residual(30) = lhs - rhs;
lhs = y(21);
rhs = y(21)*params(13)-x(1);
residual(31) = lhs - rhs;
lhs = y(32);
rhs = y(22);
residual(32) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
