function [y, T] = static_3(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(12)=100*log(y(4));
  y(13)=100*log(y(7));
  y(14)=100*log(y(8));
  y(11)=100*log(y(2));
  y(10)=100*log(y(3));
  y(9)=100*log(y(1));
  y(6)=y(1)*params(6)*4/y(3);
end
