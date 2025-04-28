function [y, T, residual, g1] = static_6(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(8, 1);
  T(1)=params(5)/(1-params(3))/(params(5)/(1-params(3))/(params(6)/(1-params(3))))^(1-params(10));
  residual(1)=(y(8))-(T(1)*y(6)^(-params(10)));
  T(2)=(1-params(3))*params(2)*exp(y(10));
  T(3)=y(3)^params(7);
  T(4)=((params(5)/(1-params(3))/(params(6)/(1-params(3)))*params(4)/(params(5)+params(4))-(1-params(3))*(params(5)/(1-params(3))/(params(6)/(1-params(3)))*params(4)/(params(5)+params(4))*(1-params(6)/(1-params(3)))+(params(4)-params(3))/(1-params(3))*(1-params(4)/(params(5)+params(4)))))*(((1-params(2))/params(2)+params(3))/((params(1)-params(9))*params(6)*(1-params(8))/(params(6)*(1-params(8))+(1-params(2))/params(2)+params(4)+params(5)/(1-params(3))*params(8))*(1+(1-params(2))/params(2))))^(1/params(7)))^params(7);
  residual(2)=(y(1))-((T(3)-T(2)*T(3))/T(4));
  T(5)=exp(y(9));
  residual(3)=(y(1)/y(8))-(T(2)*(T(5)-y(2)-y(1)+y(1)*(1-(params(4)-params(3))/(1-params(3)))/y(8)));
  residual(4)=(y(7))-(T(1)*y(6)^(1-params(10)));
  residual(5)=(y(6))-(y(4)/y(5));
  residual(6)=(y(5))-(y(5)*(1-(1-params(3))*y(7))+params(4)*(1-y(5)));
  residual(7)=(y(4))-(y(3)+(1-params(3))*(y(4)*(1-y(8))+(params(4)-params(3))/(1-params(3))*(1-y(5))));
  T(6)=params(8)*exp(y(11));
  residual(8)=(y(2))-(T(6)*(T(5)-y(1)+y(1)*y(6)/(1-params(3)))+params(9)*(1-T(6)));
  T(7)=getPowerDeriv(y(3),params(7),1);
if nargout > 3
    g1_v = NaN(21, 1);
g1_v(1)=(-(T(1)*getPowerDeriv(y(6),(-params(10)),1)));
g1_v(2)=(-(T(1)*getPowerDeriv(y(6),1-params(10),1)));
g1_v(3)=1;
g1_v(4)=(-(T(6)*y(1)*1/(1-params(3))));
g1_v(5)=(-((T(7)-T(2)*T(7))/T(4)));
g1_v(6)=(-1);
g1_v(7)=T(2);
g1_v(8)=1;
g1_v(9)=1;
g1_v(10)=(-(y(5)*(-(1-params(3)))));
g1_v(11)=(-(1/y(5)));
g1_v(12)=1-(1-params(3))*(1-y(8));
g1_v(13)=(-((-y(4))/(y(5)*y(5))));
g1_v(14)=1-(1-(1-params(3))*y(7)-params(4));
g1_v(15)=(-((1-params(3))*(-((params(4)-params(3))/(1-params(3))))));
g1_v(16)=1;
g1_v(17)=(-y(1))/(y(8)*y(8))-T(2)*(-(y(1)*(1-(params(4)-params(3))/(1-params(3)))))/(y(8)*y(8));
g1_v(18)=(-((1-params(3))*(-y(4))));
g1_v(19)=1;
g1_v(20)=1/y(8)-T(2)*((1-(params(4)-params(3))/(1-params(3)))/y(8)-1);
g1_v(21)=(-(T(6)*(y(6)/(1-params(3))-1)));
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 8, 8);
end
end
