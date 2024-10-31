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
    T = flow_specification.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(32, 1);
lhs = T(8)/y(13);
rhs = T(1)*T(11);
residual(1) = lhs - rhs;
lhs = y(17);
rhs = T(12)/T(5);
residual(2) = lhs - rhs;
lhs = y(8);
rhs = T(9)*(y(17)+y(11)*T(8)/(1-params(3)))+params(9)*params(1)*(1-T(9));
residual(3) = lhs - rhs;
lhs = y(11);
rhs = y(9)/y(10);
residual(4) = lhs - rhs;
lhs = y(12);
rhs = y(11)^(1-params(10))*T(4);
residual(5) = lhs - rhs;
lhs = y(13);
rhs = y(11)^(-params(10))*T(4);
residual(6) = lhs - rhs;
lhs = y(24);
rhs = 1-y(10);
residual(7) = lhs - rhs;
lhs = y(24);
rhs = y(26)+y(25);
residual(8) = lhs - rhs;
lhs = y(16);
rhs = y(14)^T(2);
residual(9) = lhs - rhs;
lhs = y(22);
rhs = y(21)^(-params(5));
residual(10) = lhs - rhs;
lhs = y(18);
rhs = y(45)*params(2)*(1-params(3))/y(22)*(y(43)+y(44));
residual(11) = lhs - rhs;
lhs = y(20);
rhs = y(26)*T(12);
residual(12) = lhs - rhs;
lhs = y(20);
rhs = y(21)+y(9)*T(8);
residual(13) = lhs - rhs;
lhs = y(15);
rhs = exp(y(28))*y(25)*T(7)/T(10);
residual(14) = lhs - rhs;
lhs = y(18);
rhs = y(16)*T(10)/T(5);
residual(15) = lhs - rhs;
lhs = y(23);
rhs = y(20)+y(18)*y(15);
residual(16) = lhs - rhs;
lhs = y(23);
rhs = y(17)*y(24)+y(14)*y(19);
residual(17) = lhs - rhs;
lhs = y(10);
rhs = (1-(1-params(3))*y(2))*y(1)+params(6)*(1-y(1));
residual(18) = lhs - rhs;
lhs = y(14);
rhs = (1-params(3))*(y(3)+y(4));
residual(19) = lhs - rhs;
lhs = y(27);
rhs = y(8)*y(24)/y(23);
residual(20) = lhs - rhs;
lhs = y(29);
rhs = y(28)-y(5);
residual(21) = lhs - rhs;
lhs = y(30);
rhs = y(29)+y(6)+y(7);
residual(22) = lhs - rhs;
lhs = y(31);
rhs = y(10)-(steady_state(3));
residual(23) = lhs - rhs;
lhs = y(32);
rhs = y(9)-(steady_state(2));
residual(24) = lhs - rhs;
lhs = y(33);
rhs = 100*log(y(10));
residual(25) = lhs - rhs;
lhs = y(34);
rhs = 100*log(y(9));
residual(26) = lhs - rhs;
lhs = y(35);
rhs = 100*log(y(11));
residual(27) = lhs - rhs;
lhs = y(36);
rhs = 100*log(y(14));
residual(28) = lhs - rhs;
lhs = y(37);
rhs = 100*log(y(15));
residual(29) = lhs - rhs;
lhs = y(38);
rhs = 100*log(y(21));
residual(30) = lhs - rhs;
lhs = y(28);
rhs = y(5)*params(13)-x(it_, 1);
residual(31) = lhs - rhs;
lhs = y(39);
rhs = y(6);
residual(32) = lhs - rhs;

end
