function [trend, error_status] = Func_0016_get_trend_range_monte_carlo(xval, ens_mean, ens_std)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function ensemble_name = Func_0015_GCM_CMIP6_ensname(model_name);
%
% get a ensemble name corresponding to GCM model name
%
%  input:
%  model_name             GCM model name (string)
%
%  output:
%  ensemble_name             CMIP6 ensemble name (string)
%
%  e-mail:kimyy308@snu.ac.kr
%
%  Updated    27-Jul-2021 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trend.val
trend.upper
trend.lower
trend.all

data.all
data.avg
data.std


switch model_name
    case {'CMCC-CM2-HR4', 'EC-Earth3-Veg', 'ACCESS-CM2', 'CMCC-ESM2'}
        ensemble_name='r1i1p1f1';
    case {'CNRM-ESM2-1', 'CNRM-CM6-1-HR'}
        ensemble_name='r1i1p1f2';
end

error_status=1;
end

