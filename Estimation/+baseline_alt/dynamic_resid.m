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
    T = baseline_alt.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(37, 1);
lhs = y(11)/y(18);
rhs = T(3)*T(4);
residual(1) = lhs - rhs;
lhs = y(22);
rhs = T(18)/T(2);
residual(2) = lhs - rhs;
lhs = y(12);
rhs = params(10)*(y(22)-y(11)+y(11)*y(16)/(1-params(3)))+params(11)*params(1)*(1-params(10));
residual(3) = lhs - rhs;
lhs = y(10);
rhs = T(19)^params(9);
residual(4) = lhs - rhs;
lhs = y(11);
rhs = y(10)-T(3)*y(47);
residual(5) = lhs - rhs;
lhs = y(16);
rhs = y(14)/y(15);
residual(6) = lhs - rhs;
lhs = y(17);
rhs = y(16)^(1-params(12))*T(6);
residual(7) = lhs - rhs;
lhs = y(18);
rhs = y(16)^(-params(12))*T(6);
residual(8) = lhs - rhs;
lhs = y(30);
rhs = 1-y(15);
residual(9) = lhs - rhs;
lhs = y(30);
rhs = y(32)+y(31);
residual(10) = lhs - rhs;
lhs = y(21);
rhs = y(19)^(1/(params(4)-1));
residual(11) = lhs - rhs;
lhs = y(28);
rhs = y(26)^(-params(5));
residual(12) = lhs - rhs;
lhs = y(23);
rhs = y(54)*params(2)*(1-params(3))/y(28)*(y(52)+y(53));
residual(13) = lhs - rhs;
lhs = y(25);
rhs = y(32)*T(18);
residual(14) = lhs - rhs;
lhs = y(25);
rhs = y(26)+y(27);
residual(15) = lhs - rhs;
lhs = y(27);
rhs = T(20)*T(21);
residual(16) = lhs - rhs;
lhs = y(20);
rhs = exp(y(34))*y(31)*T(14)/T(15);
residual(17) = lhs - rhs;
lhs = y(23);
rhs = y(21)*T(15)/T(2);
residual(18) = lhs - rhs;
lhs = y(29);
rhs = y(25)+y(23)*y(20);
residual(19) = lhs - rhs;
lhs = y(29);
rhs = y(22)*y(30)+y(19)*y(24);
residual(20) = lhs - rhs;
lhs = y(14);
rhs = y(13)+(1-params(3))*((1-y(4))*y(1)+(params(6)-params(3))/(1-params(3))*(1-y(2)));
residual(21) = lhs - rhs;
lhs = y(15);
rhs = y(2)*(1-(1-params(3))*y(3))+params(6)*(1-y(2));
residual(22) = lhs - rhs;
lhs = y(19);
rhs = (1-params(3))*(y(5)+y(6));
residual(23) = lhs - rhs;
lhs = y(33);
rhs = y(12)*y(30)/y(29);
residual(24) = lhs - rhs;
lhs = y(35);
rhs = y(34)-y(7);
residual(25) = lhs - rhs;
lhs = y(36);
rhs = y(35)+y(8)+y(9);
residual(26) = lhs - rhs;
lhs = y(37);
rhs = y(15)-(steady_state(6));
residual(27) = lhs - rhs;
lhs = y(38);
rhs = y(14)-(steady_state(5));
residual(28) = lhs - rhs;
lhs = y(39);
rhs = 100*log(y(15));
residual(29) = lhs - rhs;
lhs = y(40);
rhs = 100*log(y(14));
residual(30) = lhs - rhs;
lhs = y(41);
rhs = 100*log(y(16));
residual(31) = lhs - rhs;
lhs = y(42);
rhs = 100*log(y(13));
residual(32) = lhs - rhs;
lhs = y(43);
rhs = 100*log(y(19));
residual(33) = lhs - rhs;
lhs = y(44);
rhs = 100*log(y(20));
residual(34) = lhs - rhs;
lhs = y(45);
rhs = 100*log(y(26));
residual(35) = lhs - rhs;
lhs = y(34);
rhs = y(7)*params(14)-x(it_, 1);
residual(36) = lhs - rhs;
lhs = y(46);
rhs = y(8);
residual(37) = lhs - rhs;

end
