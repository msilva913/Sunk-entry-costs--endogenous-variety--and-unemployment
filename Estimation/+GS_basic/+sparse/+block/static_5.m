function [y, T] = static_5(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(22)=y(12);
  y(23)=y(14);
  y(13)=y(12)+y(12)+y(22);
  y(15)=y(14)+y(14)+y(23);
end
