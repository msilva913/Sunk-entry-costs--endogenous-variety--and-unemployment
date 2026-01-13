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
    T = BGM_translog.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(21, 1);
    residual(1) = (y(23)) - (1+1/(params(9)*y(6)));
    residual(2) = (y(24)) - (y(23)/(y(23)-1));
    residual(3) = (y(15)) - (sqrt(y(6))*T(1));
    residual(4) = (T(2)) - (T(3)*(y(27)+y(26))/y(8));
    residual(5) = (y(4)) - (y(15)*y(21)*y(11));
    residual(6) = (y(9)) - (y(15)*y(21)/y(23));
    residual(7) = (y(10)) - ((T(2)*y(9)/params(4))^params(5));
    residual(8) = (y(12)) - (y(10)-y(11));
    residual(9) = (y(8)) - (y(9)*params(6)/y(21));
    residual(10) = (y(7)) - ((y(23)-1)/y(23)*y(4)/y(6));
    residual(11) = (y(4)+y(8)*y(5)) - (y(9)*y(10)+y(6)*y(7));
    residual(12) = (y(13)) - (y(4)+y(8)*y(5));
    residual(13) = (y(14)) - (y(13)/y(15));
    residual(14) = (y(6)) - ((1-params(3))*(y(2)+y(1)));
    residual(15) = (log(y(21))) - (params(8)*log(y(3))+x(it_, 1));
    residual(16) = (y(16)) - (100*log(y(13)));
    residual(17) = (y(17)) - (100*log(y(14)));
    residual(18) = (y(18)) - (100*log(y(10)));
    residual(19) = (y(19)) - (100*log(y(11)));
    residual(20) = (y(20)) - (100*log(y(12)));
    residual(21) = (y(22)) - (log(y(21))*100);

end
