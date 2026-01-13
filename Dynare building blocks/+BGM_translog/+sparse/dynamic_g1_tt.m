function [T_order, T] = dynamic_g1_tt(y, x, params, steady_state, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = BGM_translog.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
T_order = 1;
if size(T, 1) < 5
    T = [T; NaN(5 - size(T, 1), 1)];
end
T(4) = getPowerDeriv(y(22),(-params(7)),1);
T(5) = getPowerDeriv(T(2)*y(27)/params(4),params(5),1);
end
