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
    T = GS_HRW.static_resid_tt(T, y, x, params);
end
residual = zeros(23, 1);
lhs = y(1)/y(8);
rhs = T(3);
residual(1) = lhs - rhs;
lhs = y(1);
rhs = (T(4)-(1-params(3))*params(2)*exp(y(10))*T(4))/T(5);
residual(2) = lhs - rhs;
lhs = y(2);
rhs = T(6)+params(9)*(1-params(8)*exp(y(11)));
residual(3) = lhs - rhs;
lhs = y(4);
rhs = y(3)+(1-params(3))*(y(4)*(1-y(8))+T(2)*(1-y(5)));
residual(4) = lhs - rhs;
lhs = y(5);
rhs = y(5)*(1-(1-params(3))*y(7))+params(4)*(1-y(5));
residual(5) = lhs - rhs;
lhs = y(6);
rhs = y(4)/y(5);
residual(6) = lhs - rhs;
lhs = y(7);
rhs = y(8)*y(6);
residual(7) = lhs - rhs;
lhs = y(8);
rhs = T(7)^((-1)/params(10));
residual(8) = lhs - rhs;
residual(9) = y(12);
lhs = y(13);
rhs = y(12)+y(12)+y(12);
residual(10) = lhs - rhs;
residual(11) = y(14);
lhs = y(15);
rhs = y(14)+y(14)+y(14);
residual(12) = lhs - rhs;
lhs = y(16);
rhs = y(5)-(y(5));
residual(13) = lhs - rhs;
lhs = y(17);
rhs = y(4)-(y(4));
residual(14) = lhs - rhs;
lhs = y(18);
rhs = 100*log(y(5)/(y(5)));
residual(15) = lhs - rhs;
lhs = y(19);
rhs = 100*log(y(4)/(y(4)));
residual(16) = lhs - rhs;
lhs = y(20);
rhs = 100*log(y(6)/(y(6)));
residual(17) = lhs - rhs;
lhs = y(21);
rhs = 100*log(y(3)/(y(3)));
residual(18) = lhs - rhs;
lhs = y(9);
rhs = y(9)*params(11)-x(1);
residual(19) = lhs - rhs;
lhs = y(10);
rhs = y(10)*params(12)+x(2);
residual(20) = lhs - rhs;
lhs = y(11);
rhs = y(11)*params(13)+x(3);
residual(21) = lhs - rhs;
lhs = y(22);
rhs = y(12);
residual(22) = lhs - rhs;
lhs = y(23);
rhs = y(14);
residual(23) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
