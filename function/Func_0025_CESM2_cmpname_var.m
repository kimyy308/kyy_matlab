function component = Func_0025_CESM2_cmpname_var(var_model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [component] = Func_0025_CESM2_cmpname_var(var_model);
%
%  get component (ocn, atm, lnd, ...) from variable name
%
%  input:
%  var_model    CESM2 variable name
%
%  output:
%  component    Component
%
%  e-mail:      kimyy308@pusan.ac.kr
%
%  Updated      11-Jan-2023 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    switch var_model
        case {'TEMP', 'SALT', 'SSH', 'UVEL', 'VVEL', 'diatChl', 'diazChl', 'spChl'...
               'dust_FLUX_IN', 'dust_REMIN', 'dustToSed', 'pH_3D', 'Fe', 'O2', ...
               'P_iron_FLUX_100m', 'POC_FLUX_100m', 'POP_FLUX_100m', ...
               'SiO2_FLUX_100m', 'TMXL_DR', 'XMXL_DR' }
            component='ocn'; %t.component
        case {'AODDUST', 'AODVIS', 'dst_a1SF', 'dst_a2SF', 'dst_a3SF', ...
                'PRECT', 'PSL', 'SST', 'dst_a1', 'dst_a2', 'dst_a3', ...
                'U', 'V', 'ATM_COARSE_DUST_FLUX_CPL', 'ATM_FINE_DUST_FLUX_CPL', 'SEAICE_DUST_FLUX_CPL'}
            component='atm';
        case {'DSTFLXT', 'DSTDEP', 'DSL', 'QSOIL', 'TWS'} %TOTVEGC, FIRE
            component='lnd';
        case {'TOTAL_DISCHARGE_TO_OCEAN_LIQ', ...
                'TOTAL_DISCHARGE_TO_OCEAN_ICE'}
            component='rof';
    end
end