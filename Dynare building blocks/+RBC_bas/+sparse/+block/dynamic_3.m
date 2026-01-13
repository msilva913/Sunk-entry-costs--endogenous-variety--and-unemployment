function [y, T] = dynamic_3(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(26)=100*log(y(18));
  y(27)=100*log(y(21));
  y(28)=100*log(y(22));
  y(25)=100*log(y(16));
  y(24)=100*log(y(17));
  y(23)=100*log(y(15));
  y(20)=y(15)*params(6)*4/y(3);
end
