function [y, T, residual, g1] = static_2(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(6, 1);
  T(1)=y(2)^(-params(4));
  residual(1)=(params(2)*y(4)^(1/params(3)))-(T(1)*y(7));
  residual(2)=(y(3))-(y(3)*(1-params(5))+y(8));
  residual(3)=(y(1))-(y(2)+y(8));
  T(2)=exp(y(5));
  T(3)=T(2)*y(3)^params(6);
  T(4)=y(4)^(1-params(6));
  residual(4)=(y(1))-(T(3)*T(4));
  residual(5)=(y(7))-(y(1)*(1-params(6))/y(4));
  T(5)=params(6)*T(2)*(y(3)/y(4))^(params(6)-1)+1-params(5);
  residual(6)=(T(1))-(T(1)*params(1)*T(5));
  T(6)=getPowerDeriv(y(2),(-params(4)),1);
  T(7)=getPowerDeriv(y(3)/y(4),params(6)-1,1);
if nargout > 3
    g1_v = NaN(17, 1);
g1_v(1)=(-T(1));
g1_v(2)=1;
g1_v(3)=(-1);
g1_v(4)=(-1);
g1_v(5)=(-(y(7)*T(6)));
g1_v(6)=(-1);
g1_v(7)=T(6)-T(5)*params(1)*T(6);
g1_v(8)=1;
g1_v(9)=1;
g1_v(10)=(-((1-params(6))/y(4)));
g1_v(11)=params(2)*getPowerDeriv(y(4),1/params(3),1);
g1_v(12)=(-(T(3)*getPowerDeriv(y(4),1-params(6),1)));
g1_v(13)=(-((-(y(1)*(1-params(6))))/(y(4)*y(4))));
g1_v(14)=(-(T(1)*params(1)*params(6)*T(2)*T(7)*(-y(3))/(y(4)*y(4))));
g1_v(15)=1-(1-params(5));
g1_v(16)=(-(T(4)*T(2)*getPowerDeriv(y(3),params(6),1)));
g1_v(17)=(-(T(1)*params(1)*params(6)*T(2)*1/y(4)*T(7)));
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 6, 6);
end
end
