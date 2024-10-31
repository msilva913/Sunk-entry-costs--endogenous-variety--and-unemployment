function [y, T] = dynamic_3(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(77)=100*log(y(51));
  y(76)=100*log(y(50));
  y(78)=100*log(y(57));
  y(79)=100*log(y(60));
  y(75)=100*log(y(44));
  y(74)=100*log(y(47));
  y(73)=100*log(y(45));
  y(72)=100*log(y(46));
  y(71)=y(45)-(steady_state(5));
  y(70)=y(46)-(steady_state(6));
  y(64)=y(43)*y(61)/y(60);
end
