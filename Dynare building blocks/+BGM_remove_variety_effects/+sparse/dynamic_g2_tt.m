function [T_order, T] = dynamic_g2_tt(y, x, params, steady_state, T_order, T)
if T_order >= 2
    return
end
[T_order, T] = BGM_remove_variety_effects.sparse.dynamic_g1_tt(y, x, params, steady_state, T_order, T);
T_order = 2;
if size(T, 1) < 8
    T = [T; NaN(8 - size(T, 1), 1)];
end
T(7) = getPowerDeriv(y(20),(-params(8)),2);
T(8) = getPowerDeriv(T(1)*y(25)/params(4),params(5),2);
end
