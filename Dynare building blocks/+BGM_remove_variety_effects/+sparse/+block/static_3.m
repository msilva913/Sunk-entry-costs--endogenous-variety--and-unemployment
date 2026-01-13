function [y, T] = static_3(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(19)=T(1)*100;
  y(6)=1.0*y(18)/(params(6)/(params(6)-1));
  y(5)=y(6)*params(7)/y(18);
end
