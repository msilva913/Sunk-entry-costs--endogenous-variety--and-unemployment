function ds = dynamic_set_auxiliary_series(ds, params)
%
% Status : Computes Auxiliary variables of the dynamic model and returns a dseries
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

ds.AUX_ENDO_LAG_25_1=ds.cons_share_obs_m(-1);
ds.AUX_ENDO_LAG_33_1=ds.x_obs_m(-1);
ds.AUX_ENDO_LAG_35_1=ds.w_obs_m(-1);
end
