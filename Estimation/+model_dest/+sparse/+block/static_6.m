function [y, T] = static_6(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(36)=100*log(y(10));
  y(37)=100*log(y(11));
  y(18)=T(5)/(1+params(9))*(y(4)/T(5))^(1+params(9));
end
