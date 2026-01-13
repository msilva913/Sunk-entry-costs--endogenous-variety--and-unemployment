function [T_order, T] = dynamic_g1_tt(y, x, params, steady_state, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = RBC_bas.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
T_order = 1;
if size(T, 1) < 8
    T = [T; NaN(8 - size(T, 1), 1)];
end
T(7) = getPowerDeriv(y(16),(-params(4)),1);
T(8) = getPowerDeriv(y(17)/y(32),params(6)-1,1);
end
