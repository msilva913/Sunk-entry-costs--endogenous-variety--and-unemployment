function [y, T] = dynamic_3(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(44)=100*log(y(26)/(steady_state(3)));
  y(37)=log(y(25)/y(2));
  y(46)=y(14);
  y(43)=100*log(y(29)/(steady_state(6)));
  y(42)=100*log(y(27)/(steady_state(4)));
  y(41)=100*log(y(28)/(steady_state(5)));
  y(40)=y(27)-(steady_state(4));
  y(39)=y(28)-(steady_state(5));
  y(38)=y(37)+y(14)+y(23);
end
