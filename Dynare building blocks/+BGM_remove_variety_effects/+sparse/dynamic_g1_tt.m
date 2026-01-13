function [T_order, T] = dynamic_g1_tt(y, x, params, steady_state, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = BGM_remove_variety_effects.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
T_order = 1;
if size(T, 1) < 6
    T = [T; NaN(6 - size(T, 1), 1)];
end
T(3) = getPowerDeriv(y(20),(-params(8)),1);
T(4) = y(25)*T(3)/params(4);
T(5) = getPowerDeriv(T(1)*y(25)/params(4),params(5),1);
T(6) = params(1)*(1-params(3))*getPowerDeriv(y(39),(-params(8)),1);
end
