function [y, T, residual, g1] = dynamic_2(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(7, 1);
  y(29)=y(27)/y(28);
  T(1)=params(5)/(1-params(3))/(params(5)/(1-params(3))/(params(6)/(1-params(3))))^(1-params(10));
  residual(1)=(y(30))-(y(29)^(1-params(10))*T(1));
  residual(2)=(y(27))-(y(26)+(1-params(3))*((1-y(8))*y(4)+(params(4)-params(3))/(1-params(3))*(1-y(5))));
  residual(3)=(y(28))-(y(5)*(1-(1-params(3))*y(7))+params(4)*(1-y(5)));
  residual(4)=(y(31))-(y(29)^(-params(10))*T(1));
  T(2)=(1-params(3))*params(2)*exp(y(33));
  residual(5)=(y(24)/y(31))-(T(2)*(exp(y(55))-y(48)-y(47)+y(47)*(1-(params(4)-params(3))/(1-params(3)))/y(54)));
  T(3)=((params(5)/(1-params(3))/(params(6)/(1-params(3)))*params(4)/(params(5)+params(4))-(1-params(3))*(params(5)/(1-params(3))/(params(6)/(1-params(3)))*params(4)/(params(5)+params(4))*(1-params(6)/(1-params(3)))+(params(4)-params(3))/(1-params(3))*(1-params(4)/(params(5)+params(4)))))*(((1-params(2))/params(2)+params(3))/((params(1)-params(9))*params(6)*(1-params(8))/(params(6)*(1-params(8))+(1-params(2))/params(2)+params(4)+params(5)/(1-params(3))*params(8))*(1+(1-params(2))/params(2))))^(1/params(7)))^params(7);
  residual(6)=(y(24))-((y(26)^params(7)-T(2)*y(49)^params(7))/T(3));
  T(4)=params(8)*exp(y(34));
  residual(7)=(y(25))-(T(4)*(exp(y(32))-y(24)+y(24)*y(29)/(1-params(3)))+params(9)*(1-T(4)));
  T(5)=getPowerDeriv(y(29),1-params(10),1);
  T(6)=getPowerDeriv(y(29),(-params(10)),1);
if nargout > 3
    g1_v = NaN(26, 1);
g1_v(1)=(-(y(5)*(-(1-params(3)))));
g1_v(2)=(-((1-params(3))*(1-y(8))));
g1_v(3)=(-((1-params(3))*(-((params(4)-params(3))/(1-params(3))))));
g1_v(4)=(-(1-(1-params(3))*y(7)-params(4)));
g1_v(5)=(-((1-params(3))*(-y(4))));
g1_v(6)=1;
g1_v(7)=(-(T(1)*1/y(28)*T(5)));
g1_v(8)=1;
g1_v(9)=(-(T(1)*1/y(28)*T(6)));
g1_v(10)=(-(T(4)*y(24)*1/y(28)/(1-params(3))));
g1_v(11)=(-(T(1)*(-y(27))/(y(28)*y(28))*T(5)));
g1_v(12)=1;
g1_v(13)=(-(T(1)*(-y(27))/(y(28)*y(28))*T(6)));
g1_v(14)=(-(T(4)*y(24)*(-y(27))/(y(28)*y(28))/(1-params(3))));
g1_v(15)=1;
g1_v(16)=(-y(24))/(y(31)*y(31));
g1_v(17)=1/y(31);
g1_v(18)=1;
g1_v(19)=(-(T(4)*(y(29)/(1-params(3))-1)));
g1_v(20)=(-1);
g1_v(21)=(-(getPowerDeriv(y(26),params(7),1)/T(3)));
g1_v(22)=1;
g1_v(23)=(-(T(2)*(-(y(47)*(1-(params(4)-params(3))/(1-params(3)))))/(y(54)*y(54))));
g1_v(24)=(-(T(2)*((1-(params(4)-params(3))/(1-params(3)))/y(54)-1)));
g1_v(25)=(-((-(T(2)*getPowerDeriv(y(49),params(7),1)))/T(3)));
g1_v(26)=T(2);
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 7, 21);
end
end
