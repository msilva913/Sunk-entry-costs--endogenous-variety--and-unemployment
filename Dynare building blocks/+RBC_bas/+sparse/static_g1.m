function [g1, T_order, T] = static_g1(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T_order, T)
if nargin < 8
    T_order = -1;
    T = NaN(7, 1);
end
[T_order, T] = RBC_bas.sparse.static_g1_tt(y, x, params, T_order, T);
g1_v = NaN(35, 1);
g1_v(1)=1;
g1_v(2)=1;
g1_v(3)=(-((1-params(6))/y(4)));
g1_v(4)=(-(params(6)*4/y(3)));
g1_v(5)=(-(100*1/y(1)));
g1_v(6)=T(6)-T(3)*params(1)*T(6);
g1_v(7)=(-(y(7)*T(6)));
g1_v(8)=(-1);
g1_v(9)=(-(100*1/y(2)));
g1_v(10)=(-(T(1)*params(1)*params(6)*exp(y(5))*1/y(4)*T(7)));
g1_v(11)=1-(1-params(5));
g1_v(12)=(-(T(5)*exp(y(5))*getPowerDeriv(y(3),params(6),1)));
g1_v(13)=(-((-(y(1)*params(6)*4))/(y(3)*y(3))));
g1_v(14)=(-(100*1/y(3)));
g1_v(15)=(-(T(1)*params(1)*params(6)*exp(y(5))*T(7)*(-y(3))/(y(4)*y(4))));
g1_v(16)=params(2)*getPowerDeriv(y(4),1/params(3),1);
g1_v(17)=(-(T(4)*getPowerDeriv(y(4),1-params(6),1)));
g1_v(18)=(-((-(y(1)*(1-params(6))))/(y(4)*y(4))));
g1_v(19)=(-(100*1/y(4)));
g1_v(20)=(-(T(1)*params(1)*T(2)));
g1_v(21)=(-(T(4)*T(5)));
g1_v(22)=1-params(7);
g1_v(23)=1;
g1_v(24)=(-T(1));
g1_v(25)=1;
g1_v(26)=(-(100*1/y(7)));
g1_v(27)=(-1);
g1_v(28)=(-1);
g1_v(29)=(-(100*1/y(8)));
g1_v(30)=1;
g1_v(31)=1;
g1_v(32)=1;
g1_v(33)=1;
g1_v(34)=1;
g1_v(35)=1;
if ~isoctave && matlab_ver_less_than('9.8')
    sparse_rowval = double(sparse_rowval);
    sparse_colval = double(sparse_colval);
end
g1 = sparse(sparse_rowval, sparse_colval, g1_v, 14, 14);
end
