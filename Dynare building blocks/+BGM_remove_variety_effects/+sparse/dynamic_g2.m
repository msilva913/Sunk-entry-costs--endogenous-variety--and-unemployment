function [g2_v, T_order, T] = dynamic_g2(y, x, params, steady_state, T_order, T)
if nargin < 6
    T_order = -1;
    T = NaN(8, 1);
end
[T_order, T] = BGM_remove_variety_effects.sparse.dynamic_g2_tt(y, x, params, steady_state, T_order, T);
g2_v = NaN(28, 1);
g2_v(1)=T(7);
g2_v(2)=(-((y(43)+y(42))*params(1)*(1-params(3))*getPowerDeriv(y(39),(-params(8)),2)/y(24)));
g2_v(3)=(-(T(6)/y(24)));
g2_v(4)=(-((-((y(43)+y(42))*T(6)))/(y(24)*y(24))));
g2_v(5)=(-(T(6)/y(24)));
g2_v(6)=(-((-T(2))/(y(24)*y(24))));
g2_v(7)=(-((-((-(T(2)*(y(43)+y(42))))*(y(24)+y(24))))/(y(24)*y(24)*y(24)*y(24))));
g2_v(8)=(-((-T(2))/(y(24)*y(24))));
g2_v(9)=(-1.0);
g2_v(10)=(-(T(5)*y(25)*T(7)/params(4)+T(4)*T(4)*T(8)));
g2_v(11)=(-(T(5)*T(3)/params(4)+T(4)*T(1)/params(4)*T(8)));
g2_v(12)=(-(T(1)/params(4)*T(1)/params(4)*T(8)));
g2_v(13)=(-((-params(7))/(y(37)*y(37))));
g2_v(14)=(-((-((-(y(25)*params(7)))*(y(37)+y(37))))/(y(37)*y(37)*y(37)*y(37))));
g2_v(15)=(-((-params(6))/(params(6)*y(22)*params(6)*y(22))));
g2_v(16)=(-((-((-(params(6)*y(20)))*(params(6)*params(6)*y(22)+params(6)*params(6)*y(22))))/(params(6)*y(22)*params(6)*y(22)*params(6)*y(22)*params(6)*y(22))));
g2_v(17)=1;
g2_v(18)=(-1);
g2_v(19)=(-1);
g2_v(20)=(-1);
g2_v(21)=(-(params(9)*(-1)/(y(18)*y(18))));
g2_v(22)=(-1)/(y(37)*y(37));
g2_v(23)=(-(100*(-1)/(y(29)*y(29))));
g2_v(24)=(-(100*(-1)/(y(30)*y(30))));
g2_v(25)=(-(100*(-1)/(y(26)*y(26))));
g2_v(26)=(-(100*(-1)/(y(27)*y(27))));
g2_v(27)=(-(100*(-1)/(y(28)*y(28))));
g2_v(28)=(-(100*(-1)/(y(37)*y(37))));
end
