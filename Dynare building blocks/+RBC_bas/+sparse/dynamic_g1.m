function [g1, T_order, T] = dynamic_g1(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T_order, T)
if nargin < 9
    T_order = -1;
    T = NaN(8, 1);
end
[T_order, T] = RBC_bas.sparse.dynamic_g1_tt(y, x, params, steady_state, T_order, T);
g1_v = NaN(39, 1);
g1_v(1)=(-(1-params(5)));
g1_v(2)=(-(T(6)*exp(y(19))*getPowerDeriv(y(3),params(6),1)));
g1_v(3)=(-((-(y(15)*params(6)*4))/(y(3)*y(3))));
g1_v(4)=(-params(7));
g1_v(5)=1;
g1_v(6)=1;
g1_v(7)=(-((1-params(6))/y(18)));
g1_v(8)=(-(params(6)*4/y(3)));
g1_v(9)=(-(100*1/y(15)));
g1_v(10)=T(7);
g1_v(11)=(-(y(21)*T(7)));
g1_v(12)=(-1);
g1_v(13)=(-(100*1/y(16)));
g1_v(14)=(-(T(2)*params(6)*exp(y(33))*1/y(32)*T(8)));
g1_v(15)=1;
g1_v(16)=(-(100*1/y(17)));
g1_v(17)=params(2)*getPowerDeriv(y(18),1/params(3),1);
g1_v(18)=(-(T(5)*getPowerDeriv(y(18),1-params(6),1)));
g1_v(19)=(-((-(y(15)*(1-params(6))))/(y(18)*y(18))));
g1_v(20)=(-(100*1/y(18)));
g1_v(21)=(-(T(5)*T(6)));
g1_v(22)=1;
g1_v(23)=1;
g1_v(24)=(-T(1));
g1_v(25)=1;
g1_v(26)=(-(100*1/y(21)));
g1_v(27)=(-1);
g1_v(28)=(-1);
g1_v(29)=(-(100*1/y(22)));
g1_v(30)=1;
g1_v(31)=1;
g1_v(32)=1;
g1_v(33)=1;
g1_v(34)=1;
g1_v(35)=1;
g1_v(36)=(-(T(4)*params(1)*getPowerDeriv(y(30),(-params(4)),1)));
g1_v(37)=(-(T(2)*params(6)*exp(y(33))*T(8)*(-y(17))/(y(32)*y(32))));
g1_v(38)=(-(T(2)*T(3)));
g1_v(39)=(-1);
if ~isoctave && matlab_ver_less_than('9.8')
    sparse_rowval = double(sparse_rowval);
    sparse_colval = double(sparse_colval);
end
g1 = sparse(sparse_rowval, sparse_colval, g1_v, 14, 43);
end
