function [T_order, T] = static_g1_tt(y, x, params, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = RBC_bas.sparse.static_resid_tt(y, x, params, T_order, T);
T_order = 1;
if size(T, 1) < 7
    T = [T; NaN(7 - size(T, 1), 1)];
end
T(6) = getPowerDeriv(y(2),(-params(4)),1);
T(7) = getPowerDeriv(y(3)/y(4),params(6)-1,1);
end
