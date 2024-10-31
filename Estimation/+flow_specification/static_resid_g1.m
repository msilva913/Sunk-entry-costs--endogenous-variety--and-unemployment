function [residual, g1] = static_resid_g1(T, y, x, params, T_flag)
% function [residual, g1] = static_resid_g1(T, y, x, params, T_flag)
%
% Wrapper function automatically created by Dynare
%

    if T_flag
        T = flow_specification.static_g1_tt(T, y, x, params);
    end
    residual = flow_specification.static_resid(T, y, x, params, false);
    g1       = flow_specification.static_g1(T, y, x, params, false);

end
