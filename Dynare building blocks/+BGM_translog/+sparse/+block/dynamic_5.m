function [y, T] = dynamic_5(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(30)=y(28)-y(29);
  y(38)=100*log(y(30));
  y(37)=100*log(y(29));
  y(31)=y(22)+y(26)*y(23);
  y(42)=y(41)/(y(41)-1);
  y(32)=y(31)/y(33);
  y(36)=100*log(y(28));
  y(35)=100*log(y(32));
  y(34)=100*log(y(31));
end
