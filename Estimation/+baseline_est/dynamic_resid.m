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
    T = baseline_est.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(49, 1);
lhs = y(22)/y(29);
rhs = T(2)*T(4);
residual(1) = lhs - rhs;
lhs = y(33);
rhs = y(32)*exp(y(48))/T(8);
residual(2) = lhs - rhs;
lhs = y(23);
rhs = T(3)*T(17)+params(10)*T(10)*(1-T(17));
residual(3) = lhs - rhs;
lhs = y(21);
rhs = T(18)^params(9);
residual(4) = lhs - rhs;
lhs = y(22);
rhs = y(21)-T(2)*y(70);
residual(5) = lhs - rhs;
lhs = y(27);
rhs = y(25)/y(26);
residual(6) = lhs - rhs;
lhs = y(28);
rhs = T(20);
residual(7) = lhs - rhs;
lhs = y(29);
rhs = T(21);
residual(8) = lhs - rhs;
lhs = y(42);
rhs = 1-y(26);
residual(9) = lhs - rhs;
lhs = y(42);
rhs = y(44)+y(43);
residual(10) = lhs - rhs;
lhs = y(32);
rhs = y(30)^T(5);
residual(11) = lhs - rhs;
lhs = y(39);
rhs = y(37)^(-params(5));
residual(12) = lhs - rhs;
lhs = y(34);
rhs = y(77)*(1-params(3))*params(2)*exp(y(49))/y(39)*(y(75)+y(76));
residual(13) = lhs - rhs;
lhs = y(36);
rhs = y(32)*exp(y(48))*y(44);
residual(14) = lhs - rhs;
lhs = y(36);
rhs = y(37)+y(38);
residual(15) = lhs - rhs;
lhs = y(38);
rhs = T(22)*T(23);
residual(16) = lhs - rhs;
lhs = y(31);
rhs = exp(y(48))*y(43)/T(16);
residual(17) = lhs - rhs;
lhs = y(34);
rhs = y(32)*T(16)/T(8);
residual(18) = lhs - rhs;
lhs = y(36)+y(34)*y(31);
rhs = y(33)*y(42)+y(30)*y(35);
residual(19) = lhs - rhs;
lhs = y(40);
rhs = y(36)+y(34)*y(31);
residual(20) = lhs - rhs;
lhs = y(41);
rhs = y(40)/y(32);
residual(21) = lhs - rhs;
lhs = y(25);
rhs = y(24)+(1-params(3))*((1-y(5))*y(2)+(params(6)-params(3))/(1-params(3))*(1-y(3)));
residual(22) = lhs - rhs;
lhs = y(26);
rhs = y(3)*(1-(1-params(3))*y(4))+params(6)*(1-y(3));
residual(23) = lhs - rhs;
lhs = y(30);
rhs = (1-params(3))*(y(6)+y(7));
residual(24) = lhs - rhs;
lhs = y(53);
rhs = y(41)/y(42);
residual(25) = lhs - rhs;
lhs = y(45);
rhs = y(23)*y(42)/y(40);
residual(26) = lhs - rhs;
lhs = y(46);
rhs = y(37)/y(40);
residual(27) = lhs - rhs;
lhs = y(47);
rhs = (y(46)+y(9)+y(18))/3;
residual(28) = lhs - rhs;
lhs = y(54);
rhs = log(y(53))-log(y(15));
residual(29) = lhs - rhs;
lhs = y(56);
rhs = log(y(23)/y(32))-log(y(1)/y(8));
residual(30) = lhs - rhs;
lhs = y(55);
rhs = y(54)+y(16)+y(19);
residual(31) = lhs - rhs;
lhs = y(57);
rhs = y(56)+y(17)+y(20);
residual(32) = lhs - rhs;
lhs = y(58);
rhs = y(26)-(steady_state(6));
residual(33) = lhs - rhs;
lhs = y(59);
rhs = y(25)-(steady_state(5));
residual(34) = lhs - rhs;
lhs = y(60);
rhs = 100*log(y(26));
residual(35) = lhs - rhs;
lhs = y(61);
rhs = 100*log(y(25));
residual(36) = lhs - rhs;
lhs = y(62);
rhs = 100*log(y(27));
residual(37) = lhs - rhs;
lhs = y(63);
rhs = 100*log(y(24));
residual(38) = lhs - rhs;
lhs = y(64);
rhs = 100*log(y(30));
residual(39) = lhs - rhs;
lhs = y(65);
rhs = 100*log(y(31));
residual(40) = lhs - rhs;
lhs = y(66);
rhs = 100*log(y(37));
residual(41) = lhs - rhs;
lhs = y(48);
rhs = params(14)*y(10)-x(it_, 1);
residual(42) = lhs - rhs;
lhs = y(49);
rhs = params(15)*y(11)+x(it_, 2);
residual(43) = lhs - rhs;
lhs = y(50);
rhs = params(16)*y(12)+x(it_, 3);
residual(44) = lhs - rhs;
lhs = y(51);
rhs = params(17)*y(13)+x(it_, 4);
residual(45) = lhs - rhs;
lhs = y(52);
rhs = params(18)*y(14)+x(it_, 5);
residual(46) = lhs - rhs;
lhs = y(67);
rhs = y(9);
residual(47) = lhs - rhs;
lhs = y(68);
rhs = y(16);
residual(48) = lhs - rhs;
lhs = y(69);
rhs = y(17);
residual(49) = lhs - rhs;

end
