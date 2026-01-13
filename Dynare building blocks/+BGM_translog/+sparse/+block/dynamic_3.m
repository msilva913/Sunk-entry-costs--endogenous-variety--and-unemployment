function [y, T, residual, g1] = dynamic_3(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(5, 1);
  T(2)=sqrt(y(24));
  T(3)=exp(1/(y(24)*2*params(9)));
  y(33)=T(2)*T(3);
  y(41)=1+1/(params(9)*y(24));
  y(27)=y(33)*y(39)/y(41);
  T(4)=y(22)^(-params(7));
  y(28)=(T(4)*y(27)/params(4))^params(5);
  residual(1)=(y(24))-((1-params(3))*(y(3)+y(2)));
  residual(2)=(y(22)+y(26)*y(23))-(y(27)*y(28)+y(24)*y(25));
  T(5)=params(1)*(1-params(3))*y(43)^(-params(7));
  T(6)=T(5)*(y(47)+y(46));
  residual(3)=(T(4))-(T(6)/y(26));
  residual(4)=(y(26))-(y(27)*params(6)/y(39));
  residual(5)=(y(25))-((y(41)-1)/y(41)*y(22)/y(24));
  T(7)=getPowerDeriv(y(22),(-params(7)),1);
  T(8)=getPowerDeriv(T(4)*y(27)/params(4),params(5),1);
  T(9)=(y(41)*y(39)*(T(3)*1/(T(2)+T(2))+T(2)*T(3)*(-(2*params(9)))/(y(24)*2*params(9)*y(24)*2*params(9)))-y(33)*y(39)*(-params(9))/(params(9)*y(24)*params(9)*y(24)))/(y(41)*y(41));
if nargout > 3
    g1_v = NaN(18, 1);
g1_v(1)=(-(1-params(3)));
g1_v(2)=(-(1-params(3)));
g1_v(3)=1;
g1_v(4)=(-(y(25)+y(28)*T(9)+y(27)*T(8)*T(4)*T(9)/params(4)));
g1_v(5)=(-(params(6)*T(9)/y(39)));
g1_v(6)=(-((y(41)-1)/y(41)*(-y(22))/(y(24)*y(24))+y(22)/y(24)*(y(41)*(-params(9))/(params(9)*y(24)*params(9)*y(24))-(y(41)-1)*(-params(9))/(params(9)*y(24)*params(9)*y(24)))/(y(41)*y(41))));
g1_v(7)=y(26);
g1_v(8)=1-y(27)*y(27)*T(7)/params(4)*T(8);
g1_v(9)=T(7);
g1_v(10)=(-((y(41)-1)/y(41)*1/y(24)));
g1_v(11)=y(23);
g1_v(12)=(-((-T(6))/(y(26)*y(26))));
g1_v(13)=1;
g1_v(14)=(-y(24));
g1_v(15)=1;
g1_v(16)=(-((y(47)+y(46))*params(1)*(1-params(3))*getPowerDeriv(y(43),(-params(7)),1)/y(26)));
g1_v(17)=(-(T(5)/y(26)));
g1_v(18)=(-(T(5)/y(26)));
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 5, 15);
end
end
