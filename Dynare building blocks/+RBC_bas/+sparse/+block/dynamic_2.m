function [y, T, residual, g1] = dynamic_2(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(6, 1);
  T(1)=y(16)^(-params(4));
  residual(1)=(params(2)*y(18)^(1/params(3)))-(T(1)*y(21));
  residual(2)=(y(15))-(y(16)+y(22));
  residual(3)=(y(21))-(y(15)*(1-params(6))/y(18));
  residual(4)=(y(17))-((1-params(5))*y(3)+y(22));
  T(2)=exp(y(19));
  T(3)=T(2)*y(3)^params(6);
  T(4)=y(18)^(1-params(6));
  residual(5)=(y(15))-(T(3)*T(4));
  T(5)=params(1)*y(30)^(-params(4));
  T(6)=params(6)*exp(y(33));
  T(7)=T(6)*(y(17)/y(32))^(params(6)-1)+1-params(5);
  residual(6)=(T(1))-(T(5)*T(7));
  T(8)=getPowerDeriv(y(16),(-params(4)),1);
  T(9)=getPowerDeriv(y(17)/y(32),params(6)-1,1);
if nargout > 3
    g1_v = NaN(19, 1);
g1_v(1)=(-(1-params(5)));
g1_v(2)=(-(T(4)*T(2)*getPowerDeriv(y(3),params(6),1)));
g1_v(3)=(-T(1));
g1_v(4)=1;
g1_v(5)=(-1);
g1_v(6)=(-1);
g1_v(7)=1;
g1_v(8)=(-((1-params(6))/y(18)));
g1_v(9)=1;
g1_v(10)=1;
g1_v(11)=(-(T(5)*T(6)*1/y(32)*T(9)));
g1_v(12)=params(2)*getPowerDeriv(y(18),1/params(3),1);
g1_v(13)=(-((-(y(15)*(1-params(6))))/(y(18)*y(18))));
g1_v(14)=(-(T(3)*getPowerDeriv(y(18),1-params(6),1)));
g1_v(15)=(-(y(21)*T(8)));
g1_v(16)=(-1);
g1_v(17)=T(8);
g1_v(18)=(-(T(5)*T(6)*T(9)*(-y(17))/(y(32)*y(32))));
g1_v(19)=(-(T(7)*params(1)*getPowerDeriv(y(30),(-params(4)),1)));
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 6, 18);
end
end
