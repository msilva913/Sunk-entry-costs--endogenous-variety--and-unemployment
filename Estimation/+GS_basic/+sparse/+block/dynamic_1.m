function [y, T] = dynamic_1(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(32)=y(9)*params(11)-x(1);
  y(33)=params(12)*y(10)+x(2);
  y(34)=params(13)*y(11)+x(3);
  y(35)=y(32)-y(9);
  y(45)=y(12);
  y(36)=y(35)+y(12)+y(22);
end
