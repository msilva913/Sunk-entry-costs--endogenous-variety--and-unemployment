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
    T = GS.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(23, 1);
lhs = y(12)/y(19);
rhs = T(2)*(exp(y(39))-y(36)-y(35)+y(35)*T(4)/y(38));
residual(1) = lhs - rhs;
lhs = y(12);
rhs = (y(14)^params(7)-T(6))/T(10);
residual(2) = lhs - rhs;
lhs = y(13);
rhs = T(3)+params(9)*(1-params(8)*exp(y(23)));
residual(3) = lhs - rhs;
lhs = y(15);
rhs = y(14)+(1-params(3)*exp(y(9)))*((1-y(4))*y(1)+T(1)*(1-y(2)));
residual(4) = lhs - rhs;
lhs = y(16);
rhs = y(2)*(1-(1-params(3)*exp(y(9)))*y(3))+(1-y(2))*y(5);
residual(5) = lhs - rhs;
lhs = y(17);
rhs = y(15)/y(16);
residual(6) = lhs - rhs;
lhs = y(18);
rhs = y(17)^(1-params(10))*T(9);
residual(7) = lhs - rhs;
lhs = y(19);
rhs = y(17)^(-params(10))*T(9);
residual(8) = lhs - rhs;
lhs = y(20);
rhs = 1-(1-params(3)*exp(y(24)))*T(4);
residual(9) = lhs - rhs;
lhs = y(25);
rhs = y(21)-y(6);
residual(10) = lhs - rhs;
lhs = y(26);
rhs = y(25)+y(10)+y(11);
residual(11) = lhs - rhs;
lhs = y(27);
rhs = y(16)-(steady_state(5));
residual(12) = lhs - rhs;
lhs = y(28);
rhs = y(15)-(steady_state(4));
residual(13) = lhs - rhs;
lhs = y(29);
rhs = y(20)-params(4);
residual(14) = lhs - rhs;
lhs = y(30);
rhs = 100*log(y(16)/(steady_state(5)));
residual(15) = lhs - rhs;
lhs = y(31);
rhs = 100*log(y(15)/(steady_state(4)));
residual(16) = lhs - rhs;
lhs = y(32);
rhs = 100*log(y(17)/(steady_state(6)));
residual(17) = lhs - rhs;
lhs = y(33);
rhs = 100*log(y(14)/(steady_state(3)));
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
lhs = y(24);
rhs = y(9)*params(14)+x(it_, 4);
residual(22) = lhs - rhs;
lhs = y(34);
rhs = y(10);
residual(23) = lhs - rhs;

end
