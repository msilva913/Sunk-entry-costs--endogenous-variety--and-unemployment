function [y, T, residual, g1] = static_3(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(8, 1);
  residual(1)=(y(3))-((1-params(3))*(y(3)+y(2)));
  residual(2)=(y(12))-(y(3)^(1/(params(6)-1)));
  residual(3)=(y(6))-(y(12)*y(18)/(params(6)/(params(6)-1)));
  T(2)=y(1)^(-params(8));
  residual(4)=(y(7))-((T(2)*y(6)/params(4))^params(5));
  residual(5)=(y(5))-(y(6)*params(7)/y(18));
  residual(6)=(y(4))-(y(1)/(params(6)*y(3)));
  residual(7)=(y(1)+y(5)*y(2))-(y(7)*y(6)+y(4)*y(3));
  residual(8)=(T(2))-(T(2)*params(1)*(1-params(3))*(y(5)+y(4))/y(5));
  T(3)=getPowerDeriv(y(1),(-params(8)),1);
  T(4)=getPowerDeriv(T(2)*y(6)/params(4),params(5),1);
if nargout > 3
    g1_v = NaN(24, 1);
g1_v(1)=(-(1-params(3)));
g1_v(2)=y(5);
g1_v(3)=1;
g1_v(4)=(-(y(18)/(params(6)/(params(6)-1))));
g1_v(5)=1;
g1_v(6)=(-(T(4)*T(2)/params(4)));
g1_v(7)=(-(params(7)/y(18)));
g1_v(8)=(-y(7));
g1_v(9)=1;
g1_v(10)=(-y(6));
g1_v(11)=1;
g1_v(12)=y(2);
g1_v(13)=(-((T(2)*params(1)*(1-params(3))*y(5)-T(2)*params(1)*(1-params(3))*(y(5)+y(4)))/(y(5)*y(5))));
g1_v(14)=1-(1-params(3));
g1_v(15)=(-(getPowerDeriv(y(3),1/(params(6)-1),1)));
g1_v(16)=(-((-(params(6)*y(1)))/(params(6)*y(3)*params(6)*y(3))));
g1_v(17)=(-y(4));
g1_v(18)=(-(y(6)*T(3)/params(4)*T(4)));
g1_v(19)=(-(1/(params(6)*y(3))));
g1_v(20)=1;
g1_v(21)=T(3)-(y(5)+y(4))*params(1)*(1-params(3))*T(3)/y(5);
g1_v(22)=1;
g1_v(23)=(-y(3));
g1_v(24)=(-(T(2)*params(1)*(1-params(3))/y(5)));
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 8, 8);
end
end
