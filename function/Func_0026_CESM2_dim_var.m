function dim = Func_0026_CESM2_dim_var(var_model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [dimension] = Func_0026_CESM2_dim_var(var_model);
%
% get dimension information from variable name
%
%  input:
%  var_model    CESM2 variable name
%
%  output:
%  dim          Dimension
%
%  e-mail:      kimyy308@pusan.ac.kr
%
%  Updated      11-Jan-2023 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch var_model
    case {'TEMP', 'SALT', 'UVEL', 'VVEL', 'diatChl', 'diazChl', 'spChl'...
           'dust_FLUX_IN', 'dust_REMIN',  'pH_3D', 'Fe', 'O2', ...
           'dst_a1', 'dst_a2', 'dst_a3', 'U', 'V'}
        dim='4d'; %[x, y, z, t]
    case {'dustToSed', 'SSH', 'AODDUST', 'AODVIS', 'dst_a1SF', 'dst_a2SF', 'dst_a3SF', ...
            'PRECT', 'PSL', 'SST',  ...
            'ATM_COARSE_DUST_FLUX_CPL', 'ATM_FINE_DUST_FLUX_CPL', 'SEAICE_DUST_FLUX_CPL', ...
           'P_iron_FLUX_100m', 'POC_FLUX_100m', 'POP_FLUX_100m', ...
           'SiO2_FLUX_100m', 'TMXL_DR', 'XMXL_DR', ...
           'DSTFLXT', 'DSTDEP', 'DSL', 'QSOIL', 'TWS',...
           'TOTAL_DISCHARGE_TO_OCEAN_LIQ', ...
            'TOTAL_DISCHARGE_TO_OCEAN_ICE'}
        dim='3d'; % [x, y, t]
end

end