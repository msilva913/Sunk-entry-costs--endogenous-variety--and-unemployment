function [y, T] = dynamic_3(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(80)=100*log(y(58));
  y(78)=100*log(y(52));
  y(81)=100*log(y(61));
  y(79)=100*log(y(53));
  y(77)=100*log(y(51));
  y(76)=100*log(y(45));
  y(75)=100*log(y(48));
  y(74)=100*log(y(46));
  y(73)=100*log(y(47));
  y(72)=y(46)-(steady_state(5));
  y(71)=y(47)-(steady_state(6));
  y(65)=y(44)*y(62)/y(61);
end
