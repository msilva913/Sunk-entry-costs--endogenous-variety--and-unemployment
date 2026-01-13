function [T_order, T] = dynamic_g2_tt(y, x, params, steady_state, T_order, T)
if T_order >= 2
    return
end
[T_order, T] = BGM.sparse.dynamic_g1_tt(y, x, params, steady_state, T_order, T);
T_order = 2;
if size(T, 1) < 9
    T = [T; NaN(9 - size(T, 1), 1)];
end
T(8) = getPowerDeriv(y(20),(-params(8)),2);
T(9) = getPowerDeriv(T(2)*y(25)/params(4),params(5),2);
end
