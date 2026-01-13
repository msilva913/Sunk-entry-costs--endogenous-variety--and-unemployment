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
    T = BGM_translog.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(21, 28);
g1(1,6)=(-((-params(9))/(params(9)*y(6)*params(9)*y(6))));
g1(1,23)=1;
g1(2,23)=(-((y(23)-1-y(23))/((y(23)-1)*(y(23)-1))));
g1(2,24)=1;
g1(3,6)=(-(T(1)*1/(sqrt(y(6))+sqrt(y(6)))+sqrt(y(6))*T(1)*(-(2*params(9)))/(y(6)*2*params(9)*y(6)*2*params(9))));
g1(3,15)=1;
g1(4,4)=T(4);
g1(4,25)=(-((y(27)+y(26))*params(1)*(1-params(3))*getPowerDeriv(y(25),(-params(7)),1)/y(8)));
g1(4,26)=(-(T(3)/y(8)));
g1(4,8)=(-((-(T(3)*(y(27)+y(26))))/(y(8)*y(8))));
g1(4,27)=(-(T(3)/y(8)));
g1(5,4)=1;
g1(5,11)=(-(y(15)*y(21)));
g1(5,15)=(-(y(21)*y(11)));
g1(5,21)=(-(y(15)*y(11)));
g1(6,9)=1;
g1(6,15)=(-(y(21)/y(23)));
g1(6,21)=(-(y(15)/y(23)));
g1(6,23)=(-((-(y(15)*y(21)))/(y(23)*y(23))));
g1(7,4)=(-(y(9)*T(4)/params(4)*T(5)));
g1(7,9)=(-(T(5)*T(2)/params(4)));
g1(7,10)=1;
g1(8,10)=(-1);
g1(8,11)=1;
g1(8,12)=1;
g1(9,8)=1;
g1(9,9)=(-(params(6)/y(21)));
g1(9,21)=(-((-(y(9)*params(6)))/(y(21)*y(21))));
g1(10,4)=(-((y(23)-1)/y(23)*1/y(6)));
g1(10,6)=(-((y(23)-1)/y(23)*(-y(4))/(y(6)*y(6))));
g1(10,7)=1;
g1(10,23)=(-(y(4)/y(6)*(y(23)-(y(23)-1))/(y(23)*y(23))));
g1(11,4)=1;
g1(11,5)=y(8);
g1(11,6)=(-y(7));
g1(11,7)=(-y(6));
g1(11,8)=y(5);
g1(11,9)=(-y(10));
g1(11,10)=(-y(9));
g1(12,4)=(-1);
g1(12,5)=(-y(8));
g1(12,8)=(-y(5));
g1(12,13)=1;
g1(13,13)=(-(1/y(15)));
g1(13,14)=1;
g1(13,15)=(-((-y(13))/(y(15)*y(15))));
g1(14,1)=(-(1-params(3)));
g1(14,2)=(-(1-params(3)));
g1(14,6)=1;
g1(15,3)=(-(params(8)*1/y(3)));
g1(15,21)=1/y(21);
g1(15,28)=(-1);
g1(16,13)=(-(100*1/y(13)));
g1(16,16)=1;
g1(17,14)=(-(100*1/y(14)));
g1(17,17)=1;
g1(18,10)=(-(100*1/y(10)));
g1(18,18)=1;
g1(19,11)=(-(100*1/y(11)));
g1(19,19)=1;
g1(20,12)=(-(100*1/y(12)));
g1(20,20)=1;
g1(21,21)=(-(100*1/y(21)));
g1(21,22)=1;

end
