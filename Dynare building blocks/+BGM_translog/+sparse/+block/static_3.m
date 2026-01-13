function [y, T, residual, g1] = static_3(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(9, 1);
  T(2)=sqrt(y(3));
  T(3)=exp(1/(y(3)*2*params(9)));
  residual(1)=(y(12))-(T(2)*T(3));
  T(4)=y(1)^(-params(7));
  residual(2)=(T(4))-(T(4)*params(1)*(1-params(3))*(y(5)+y(4))/y(5));
  residual(3)=(y(3))-((1-params(3))*(y(3)+y(2)));
  residual(4)=(y(6))-(y(12)*y(18)/y(20));
  residual(5)=(y(7))-((T(4)*y(6)/params(4))^params(5));
  residual(6)=(y(20))-(1+1/(params(9)*y(3)));
  residual(7)=(y(5))-(y(6)*params(6)/y(18));
  residual(8)=(y(4))-((y(20)-1)/y(20)*y(1)/y(3));
  residual(9)=(y(1)+y(5)*y(2))-(y(6)*y(7)+y(3)*y(4));
  T(5)=getPowerDeriv(y(1),(-params(7)),1);
  T(6)=getPowerDeriv(T(4)*y(6)/params(4),params(5),1);
if nargout > 3
    g1_v = NaN(28, 1);
g1_v(1)=1;
g1_v(2)=(-(y(18)/y(20)));
g1_v(3)=(-((T(4)*params(1)*(1-params(3))*y(5)-T(4)*params(1)*(1-params(3))*(y(5)+y(4)))/(y(5)*y(5))));
g1_v(4)=1;
g1_v(5)=y(2);
g1_v(6)=(-(1-params(3)));
g1_v(7)=y(5);
g1_v(8)=(-((-(y(12)*y(18)))/(y(20)*y(20))));
g1_v(9)=1;
g1_v(10)=(-(y(1)/y(3)*(y(20)-(y(20)-1))/(y(20)*y(20))));
g1_v(11)=T(5)-(y(5)+y(4))*params(1)*(1-params(3))*T(5)/y(5);
g1_v(12)=(-(y(6)*T(5)/params(4)*T(6)));
g1_v(13)=(-((y(20)-1)/y(20)*1/y(3)));
g1_v(14)=1;
g1_v(15)=(-(T(3)*1/(T(2)+T(2))+T(2)*T(3)*(-(2*params(9)))/(y(3)*2*params(9)*y(3)*2*params(9))));
g1_v(16)=1-(1-params(3));
g1_v(17)=(-((-params(9))/(params(9)*y(3)*params(9)*y(3))));
g1_v(18)=(-((y(20)-1)/y(20)*(-y(1))/(y(3)*y(3))));
g1_v(19)=(-y(4));
g1_v(20)=1;
g1_v(21)=(-(T(6)*T(4)/params(4)));
g1_v(22)=(-(params(6)/y(18)));
g1_v(23)=(-y(7));
g1_v(24)=(-(T(4)*params(1)*(1-params(3))/y(5)));
g1_v(25)=1;
g1_v(26)=(-y(3));
g1_v(27)=1;
g1_v(28)=(-y(6));
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 9, 9);
end
end
