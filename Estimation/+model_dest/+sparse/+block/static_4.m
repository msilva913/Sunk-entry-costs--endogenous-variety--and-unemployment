function [y, T] = static_4(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(41)=y(28);
  T(1)=1-params(2)*exp(y(27));
  y(25)=1-T(1)*(1-(params(5)-params(2))/(1-params(2)));
  y(29)=y(28)+y(28)+y(41);
end
