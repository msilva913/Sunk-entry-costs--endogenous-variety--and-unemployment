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
    residual(1) = (y(13)/y(20)) - (T(3));
    residual(2) = (y(13)) - ((y(15)^params(7)-T(4))/T(8));
    residual(3) = (y(14)) - (T(2)+params(9)*(1-params(8)*exp(y(23))));
    residual(4) = (y(16)) - (y(15)+(1-params(3))*((1-y(5))*y(2)+T(1)*(1-y(3))));
    residual(5) = (y(17)) - (y(3)*(1-(1-params(3))*y(4))+params(4)*(1-y(3)));
    residual(6) = (y(18)) - (y(16)/y(17));
    residual(7) = (y(19)) - (y(18)^(1-params(10))*T(7));
    residual(8) = (y(20)) - (y(18)^(-params(10))*T(7));
    residual(9) = (y(24)) - (y(21)-y(6));
    residual(10) = (y(25)) - (y(24)+y(9)+y(11));
    residual(11) = (y(26)) - (log(y(14)/y(1)));
    residual(12) = (y(27)) - (y(26)+y(10)+y(12));
    residual(13) = (y(28)) - (y(17)-(steady_state(5)));
    residual(14) = (y(29)) - (y(16)-(steady_state(4)));
    residual(15) = (y(30)) - (100*log(y(17)/(steady_state(5))));
    residual(16) = (y(31)) - (100*log(y(16)/(steady_state(4))));
    residual(17) = (y(32)) - (100*log(y(18)/(steady_state(6))));
    residual(18) = (y(33)) - (100*log(y(15)/(steady_state(3))));
    residual(19) = (y(21)) - (y(6)*params(11)-x(it_, 1));
    residual(20) = (y(22)) - (params(12)*y(7)+x(it_, 2));
    residual(21) = (y(23)) - (params(13)*y(8)+x(it_, 3));
    residual(22) = (y(34)) - (y(9));
    residual(23) = (y(35)) - (y(10));

end
