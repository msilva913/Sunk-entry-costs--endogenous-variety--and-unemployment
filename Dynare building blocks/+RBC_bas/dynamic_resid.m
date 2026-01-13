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
    T = RBC_bas.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(14, 1);
    residual(1) = (T(1)) - (T(2)*T(4));
    residual(2) = (params(2)*y(6)^(1/params(3))) - (T(1)*y(9));
    residual(3) = (y(5)) - ((1-params(5))*y(1)+y(10));
    residual(4) = (y(3)) - (y(4)+y(10));
    residual(5) = (y(3)) - (T(5)*T(6));
    residual(6) = (y(9)) - (y(3)*(1-params(6))/y(6));
    residual(7) = (y(8)) - (y(3)*params(6)*4/y(1));
    residual(8) = (y(7)) - (params(7)*y(2)+x(it_, 1));
    residual(9) = (y(11)) - (100*log(y(3)));
    residual(10) = (y(12)) - (100*log(y(5)));
    residual(11) = (y(13)) - (100*log(y(4)));
    residual(12) = (y(14)) - (100*log(y(6)));
    residual(13) = (y(15)) - (100*log(y(9)));
    residual(14) = (y(16)) - (100*log(y(10)));

end
