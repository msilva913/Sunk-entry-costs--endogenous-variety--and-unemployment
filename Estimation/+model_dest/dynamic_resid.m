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
    T = model_dest.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(40, 1);
    residual(1) = (y(12)/y(19)+T(13)) - (T(3)*T(17));
    residual(2) = (y(23)) - (T(18)/T(9));
    residual(3) = (y(13)) - (T(14)*(y(23)-y(12)+T(4)*(y(12)+y(19)*T(13)))+params(14)*params(10)*(1-T(14)));
    residual(4) = (y(11)) - (T(19)^params(9));
    residual(5) = (y(12)) - (y(11)-T(3)*y(51));
    residual(6) = (y(17)) - (y(15)/y(16));
    residual(7) = (y(18)) - (y(17)^(1-params(11))*T(8));
    residual(8) = (y(19)) - (y(17)^(-params(11))*T(8));
    residual(9) = (y(31)) - (1-y(16));
    residual(10) = (y(31)) - (y(33)+y(32));
    residual(11) = (y(22)) - (y(20)^T(6));
    residual(12) = (y(29)) - (y(27)^(-params(4)));
    residual(13) = (y(24)) - (y(58)*params(1)*T(2)/y(29)*(y(56)+y(57)));
    residual(14) = (y(26)) - (y(33)*T(18));
    residual(15) = (y(26)) - (y(27)+y(28)+y(19)*y(15)*T(13));
    residual(16) = (y(28)) - (T(15)/(1+params(9))*T(19)^(1+params(9)));
    residual(17) = (y(21)) - (exp(y(36))*y(32)*T(11)/T(16));
    residual(18) = (y(24)) - (y(22)*T(16)/T(9));
    residual(19) = (y(30)) - (y(26)+y(24)*y(21));
    residual(20) = (y(30)) - (y(23)*y(31)+y(20)*y(25));
    residual(21) = (y(15)) - (y(14)+T(2)*((1-y(4))*y(1)+(params(5)-params(2))/(1-params(2))*(1-y(2))));
    residual(22) = (y(16)) - (y(2)*(1-T(2)*y(3))+(1-y(2))*y(35));
    residual(23) = (y(20)) - (T(2)*(y(5)+y(6)));
    residual(24) = (y(34)) - (y(13)*y(31)/y(30));
    residual(25) = (y(35)) - (1-T(2)*T(5));
    residual(26) = (y(38)) - (y(36)-y(7));
    residual(27) = (y(39)) - (y(38)+y(9)+y(10));
    residual(28) = (y(40)) - (y(16)-(steady_state(6)));
    residual(29) = (y(41)) - (y(15)-(steady_state(5)));
    residual(30) = (y(42)) - (100*log(y(16)));
    residual(31) = (y(43)) - (100*log(y(15)));
    residual(32) = (y(44)) - (100*log(y(17)));
    residual(33) = (y(45)) - (100*log(y(14)));
    residual(34) = (y(46)) - (100*log(y(20)));
    residual(35) = (y(47)) - (100*log(y(21)));
    residual(36) = (y(48)) - (100*log(y(27)));
    residual(37) = (y(49)) - (100*log(y(30)));
    residual(38) = (y(36)) - (y(7)*params(15)-x(it_, 1));
    residual(39) = (y(37)) - (params(16)*y(8)+x(it_, 2));
    residual(40) = (y(50)) - (y(9));

end
