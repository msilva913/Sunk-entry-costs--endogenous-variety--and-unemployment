function [y, T, residual, g1] = dynamic_4(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(4, 1);
  T(2)=y(20)^(-params(8));
  y(26)=(T(2)*y(25)/params(4))^params(5);
  residual(1)=(y(22))-((1-params(3))*(y(3)+y(2)));
  residual(2)=(y(20)+y(24)*y(21))-(y(26)*y(25)+y(23)*y(22));
  residual(3)=(y(23))-(y(20)/(params(6)*y(22)));
  T(3)=params(1)*(1-params(3))*y(39)^(-params(8));
  residual(4)=(T(2))-(T(3)*(y(43)+y(42))/y(24));
  T(4)=getPowerDeriv(y(20),(-params(8)),1);
if nargout > 3
    g1_v = NaN(13, 1);
g1_v(1)=(-(1-params(3)));
g1_v(2)=(-(1-params(3)));
g1_v(3)=1;
g1_v(4)=(-y(23));
g1_v(5)=(-((-(params(6)*y(20)))/(params(6)*y(22)*params(6)*y(22))));
g1_v(6)=y(24);
g1_v(7)=(-y(22));
g1_v(8)=1;
g1_v(9)=1-y(25)*y(25)*T(4)/params(4)*getPowerDeriv(T(2)*y(25)/params(4),params(5),1);
g1_v(10)=(-(1/(params(6)*y(22))));
g1_v(11)=T(4);
g1_v(12)=(-(T(3)/y(24)));
g1_v(13)=(-((y(43)+y(42))*params(1)*(1-params(3))*getPowerDeriv(y(39),(-params(8)),1)/y(24)));
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 4, 12);
end
end
