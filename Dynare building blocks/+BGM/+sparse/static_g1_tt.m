function [T_order, T] = static_g1_tt(y, x, params, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = BGM.sparse.static_resid_tt(y, x, params, T_order, T);
T_order = 1;
if size(T, 1) < 3
    T = [T; NaN(3 - size(T, 1), 1)];
end
T(2) = getPowerDeriv(y(1),(-params(8)),1);
T(3) = getPowerDeriv(T(1)*y(6)/params(4),params(5),1);
end
