function [T_order, T] = static_g1_tt(y, x, params, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = BGM_translog.sparse.static_resid_tt(y, x, params, T_order, T);
T_order = 1;
if size(T, 1) < 4
    T = [T; NaN(4 - size(T, 1), 1)];
end
T(3) = getPowerDeriv(y(1),(-params(7)),1);
T(4) = getPowerDeriv(T(2)*y(6)/params(4),params(5),1);
end
