function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
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
%   g1
%

if T_flag
    T = RBC_bas.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(14, 20);
g1(1,4)=T(7);
g1(1,17)=(-(T(4)*params(1)*getPowerDeriv(y(17),(-params(4)),1)));
g1(1,5)=(-(T(2)*params(6)*exp(y(19))*1/y(18)*T(8)));
g1(1,18)=(-(T(2)*params(6)*exp(y(19))*T(8)*(-y(5))/(y(18)*y(18))));
g1(1,19)=(-(T(2)*T(3)));
g1(2,4)=(-(y(9)*T(7)));
g1(2,6)=params(2)*getPowerDeriv(y(6),1/params(3),1);
g1(2,9)=(-T(1));
g1(3,1)=(-(1-params(5)));
g1(3,5)=1;
g1(3,10)=(-1);
g1(4,3)=1;
g1(4,4)=(-1);
g1(4,10)=(-1);
g1(5,3)=1;
g1(5,1)=(-(T(6)*exp(y(7))*getPowerDeriv(y(1),params(6),1)));
g1(5,6)=(-(T(5)*getPowerDeriv(y(6),1-params(6),1)));
g1(5,7)=(-(T(5)*T(6)));
g1(6,3)=(-((1-params(6))/y(6)));
g1(6,6)=(-((-(y(3)*(1-params(6))))/(y(6)*y(6))));
g1(6,9)=1;
g1(7,3)=(-(params(6)*4/y(1)));
g1(7,1)=(-((-(y(3)*params(6)*4))/(y(1)*y(1))));
g1(7,8)=1;
g1(8,2)=(-params(7));
g1(8,7)=1;
g1(8,20)=(-1);
g1(9,3)=(-(100*1/y(3)));
g1(9,11)=1;
g1(10,5)=(-(100*1/y(5)));
g1(10,12)=1;
g1(11,4)=(-(100*1/y(4)));
g1(11,13)=1;
g1(12,6)=(-(100*1/y(6)));
g1(12,14)=1;
g1(13,9)=(-(100*1/y(9)));
g1(13,15)=1;
g1(14,10)=(-(100*1/y(10)));
g1(14,16)=1;

end
