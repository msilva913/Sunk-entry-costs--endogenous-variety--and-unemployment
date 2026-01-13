function [y, T] = dynamic_3(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(38)=T(1)*100;
  y(25)=1.0*y(37)/(params(6)/(params(6)-1));
  y(24)=y(25)*params(7)/y(37);
end
