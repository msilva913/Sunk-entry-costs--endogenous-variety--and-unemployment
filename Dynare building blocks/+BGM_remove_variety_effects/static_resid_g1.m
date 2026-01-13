function [residual, g1] = static_resid_g1(T, y, x, params, T_flag)
% function [residual, g1] = static_resid_g1(T, y, x, params, T_flag)
%
% Wrapper function automatically created by Dynare
%

    if T_flag
        T = BGM_remove_variety_effects.static_g1_tt(T, y, x, params);
    end
    residual = BGM_remove_variety_effects.static_resid(T, y, x, params, false);
    g1       = BGM_remove_variety_effects.static_g1(T, y, x, params, false);

end
