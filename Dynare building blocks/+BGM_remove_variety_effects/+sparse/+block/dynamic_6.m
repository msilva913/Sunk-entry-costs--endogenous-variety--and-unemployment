function [y, T] = dynamic_6(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(35)=100*log(y(27));
  y(29)=y(20)+y(24)*y(21);
  y(28)=y(26)-y(27);
  y(36)=100*log(y(28));
  y(30)=y(29)/1.0;
  y(34)=100*log(y(26));
  y(33)=100*log(y(30));
  y(32)=100*log(y(29));
end
