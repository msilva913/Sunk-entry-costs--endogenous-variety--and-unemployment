function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = GS_basic.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(23, 1);
lhs = y(13)/y(20);
rhs = T(9);
residual(1) = lhs - rhs;
lhs = y(13);
rhs = (y(15)^params(7)-T(10))/T(14);
residual(2) = lhs - rhs;
lhs = y(14);
rhs = T(3)+params(9)*(1-params(8)*exp(y(23)));
residual(3) = lhs - rhs;
lhs = y(16);
rhs = y(15)+(1-params(3))*((1-y(5))*y(2)+T(1)*(1-y(3)));
residual(4) = lhs - rhs;
lhs = y(17);
rhs = y(3)*(1-(1-params(3))*y(4))+params(4)*(1-y(3));
residual(5) = lhs - rhs;
lhs = y(18);
rhs = y(16)/y(17);
residual(6) = lhs - rhs;
lhs = y(19);
rhs = y(18)^(1-params(10))*T(13);
residual(7) = lhs - rhs;
lhs = y(20);
rhs = y(18)^(-params(10))*T(13);
residual(8) = lhs - rhs;
lhs = y(24);
rhs = y(21)-y(6);
residual(9) = lhs - rhs;
lhs = y(25);
rhs = y(24)+y(9)+y(11);
residual(10) = lhs - rhs;
lhs = y(26);
rhs = log(T(4));
residual(11) = lhs - rhs;
lhs = y(27);
rhs = y(26)+y(10)+y(12);
residual(12) = lhs - rhs;
lhs = y(28);
rhs = y(17)-(steady_state(5));
residual(13) = lhs - rhs;
lhs = y(29);
rhs = y(16)-(steady_state(4));
residual(14) = lhs - rhs;
lhs = y(30);
rhs = 100*log(T(5));
residual(15) = lhs - rhs;
lhs = y(31);
rhs = 100*log(T(6));
residual(16) = lhs - rhs;
lhs = y(32);
rhs = 100*log(T(7));
residual(17) = lhs - rhs;
lhs = y(33);
rhs = 100*log(T(8));
residual(18) = lhs - rhs;
lhs = y(21);
rhs = y(6)*params(11)-x(it_, 1);
residual(19) = lhs - rhs;
lhs = y(22);
rhs = params(12)*y(7)+x(it_, 2);
residual(20) = lhs - rhs;
lhs = y(23);
rhs = params(13)*y(8)+x(it_, 3);
residual(21) = lhs - rhs;
lhs = y(34);
rhs = y(9);
residual(22) = lhs - rhs;
lhs = y(35);
rhs = y(10);
residual(23) = lhs - rhs;

end
