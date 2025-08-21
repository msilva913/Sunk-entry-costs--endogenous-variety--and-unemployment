function [y, T] = dynamic_1(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(67)=y(26)*params(15)-x(1);
  y(68)=params(16)*y(27)+x(2);
  T(1)=params(2)*exp(y(68));
  y(66)=1-(1-T(1))*(1-(params(5)-params(2))/(1-params(2)));
  y(69)=y(67)-y(26);
  y(82)=y(28);
  y(70)=y(69)+y(28)+y(41);
end
