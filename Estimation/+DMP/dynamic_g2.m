function g2 = dynamic_g2(T, y, x, params, steady_state, it_, T_flag)
% function g2 = dynamic_g2(T, y, x, params, steady_state, it_, T_flag)
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
%   g2
%

if T_flag
    T = DMP.dynamic_g2_tt(T, y, x, params, steady_state, it_);
end
g2_i = zeros(14,1);
g2_j = zeros(14,1);
g2_v = zeros(14,1);

g2_i(1)=1;
g2_i(2)=1;
g2_i(3)=1;
g2_i(4)=2;
g2_i(5)=3;
g2_i(6)=3;
g2_i(7)=4;
g2_i(8)=4;
g2_i(9)=4;
g2_i(10)=5;
g2_i(11)=6;
g2_i(12)=7;
g2_i(13)=8;
g2_i(14)=9;
g2_j(1)=145;
g2_j(2)=253;
g2_j(3)=271;
g2_j(4)=163;
g2_j(5)=2;
g2_j(6)=18;
g2_j(7)=74;
g2_j(8)=90;
g2_j(9)=91;
g2_j(10)=109;
g2_j(11)=109;
g2_j(12)=91;
g2_j(13)=73;
g2_j(14)=109;
g2_v(1)=(-((-T(5))*(y(9)+y(9))))/(y(9)*y(9)*y(9)*y(9));
g2_v(2)=(-(params(2)*(-((-((1-params(4))*T(5)))*(y(15)+y(15))))/(y(15)*y(15)*y(15)*y(15))));
g2_v(3)=(-(params(2)*exp(y(16))));
g2_v(4)=(-(params(7)*exp(y(10))));
g2_v(5)=1;
g2_v(6)=g2_v(5);
g2_v(7)=(-((-1)/(y(6)*y(6))));
g2_v(8)=g2_v(7);
g2_v(9)=(-((-((-y(5))*(y(6)+y(6))))/(y(6)*y(6)*y(6)*y(6))));
g2_v(10)=(-(T(4)*getPowerDeriv(y(7),1-params(9),2)));
g2_v(11)=(-(T(4)*getPowerDeriv(y(7),(-params(9)),2)));
g2_v(12)=(-(100*(-(1/(steady_state(3))*1/(steady_state(3))))/(T(1)*T(1))));
g2_v(13)=(-(100*(-(1/(steady_state(2))*1/(steady_state(2))))/(T(2)*T(2))));
g2_v(14)=(-(100*(-(1/(steady_state(4))*1/(steady_state(4))))/(T(3)*T(3))));
g2 = sparse(g2_i,g2_j,g2_v,10,289);
end
