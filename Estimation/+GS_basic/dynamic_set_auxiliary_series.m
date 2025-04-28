function ds = dynamic_set_auxiliary_series(ds, params)
%
% Computes auxiliary variables of the dynamic model
%
ds.AUX_ENDO_LAG_11_1=ds.x_obs_m(-1);
ds.AUX_ENDO_LAG_13_1=ds.w_obs_m(-1);
end
