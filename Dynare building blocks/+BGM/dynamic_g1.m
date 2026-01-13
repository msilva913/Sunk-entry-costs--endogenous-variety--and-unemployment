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
    T = BGM.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(19, 26);
g1(1,4)=T(4);
g1(1,23)=(-((y(25)+y(24))*T(7)/y(8)));
g1(1,24)=(-(T(3)/y(8)));
g1(1,8)=(-((-(T(3)*(y(25)+y(24))))/(y(8)*y(8))));
g1(1,25)=(-(T(3)/y(8)));
g1(2,6)=(-(getPowerDeriv(y(6),1/(params(6)-1),1)));
g1(2,15)=1;
g1(3,4)=1;
g1(3,11)=(-(y(15)*y(21)));
g1(3,15)=(-(y(21)*y(11)));
g1(3,21)=(-(y(15)*y(11)));
g1(4,10)=(-1);
g1(4,11)=1;
g1(4,12)=1;
g1(5,9)=1;
g1(5,15)=(-(y(21)/T(1)));
g1(5,21)=(-(y(15)/T(1)));
g1(6,4)=(-(T(5)*T(6)));
g1(6,9)=(-(T(6)*T(2)/params(4)));
g1(6,10)=1;
g1(7,8)=1;
g1(7,9)=(-(params(7)/y(21)));
g1(7,21)=(-((-(y(9)*params(7)))/(y(21)*y(21))));
g1(8,4)=(-(1/(params(6)*y(6))));
g1(8,6)=(-((-(params(6)*y(4)))/(params(6)*y(6)*params(6)*y(6))));
g1(8,7)=1;
g1(9,4)=1;
g1(9,5)=y(8);
g1(9,6)=(-y(7));
g1(9,7)=(-y(6));
g1(9,8)=y(5);
g1(9,9)=(-y(10));
g1(9,10)=(-y(9));
g1(10,4)=(-1);
g1(10,5)=(-y(8));
g1(10,8)=(-y(5));
g1(10,13)=1;
g1(11,13)=(-(1/y(15)));
g1(11,14)=1;
g1(11,15)=(-((-y(13))/(y(15)*y(15))));
g1(12,1)=(-(1-params(3)));
g1(12,2)=(-(1-params(3)));
g1(12,6)=1;
g1(13,3)=(-(params(9)*1/y(3)));
g1(13,21)=1/y(21);
g1(13,26)=(-1);
g1(14,13)=(-(100*1/y(13)));
g1(14,16)=1;
g1(15,14)=(-(100*1/y(14)));
g1(15,17)=1;
g1(16,10)=(-(100*1/y(10)));
g1(16,18)=1;
g1(17,11)=(-(100*1/y(11)));
g1(17,19)=1;
g1(18,12)=(-(100*1/y(12)));
g1(18,20)=1;
g1(19,21)=(-(100*1/y(21)));
g1(19,22)=1;

end
