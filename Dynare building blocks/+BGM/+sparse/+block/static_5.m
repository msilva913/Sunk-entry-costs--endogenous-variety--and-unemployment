function [y, T] = static_5(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(16)=100*log(y(8));
  y(10)=y(1)+y(5)*y(2);
  y(9)=y(7)-y(8);
  y(17)=100*log(y(9));
  y(11)=y(10)/y(12);
  y(15)=100*log(y(7));
  y(14)=100*log(y(11));
  y(13)=100*log(y(10));
end
