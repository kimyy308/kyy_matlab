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
    case {'TEMP', 'SALT', 'UVEL', 'VVEL', 'diatChl', 'diazChl', 'spChl', 'PD', ...
           'dust_FLUX_IN', 'dust_REMIN',  'pH_3D', 'Fe', 'O2', ...
           'dst_a1', 'dst_a2', 'dst_a3', 'U', 'V', 'ALK',  'DIC',  'NO3', 'PO4', 'SiO3', 'WVEL', ...
           'UISOP', 'VISOP', 'UE_PO4', 'VN_PO4', 'WT_PO4', ...
           'CLOUD','OMEGA', 'Q', 'RELHUM', 'diatC', 'diazC', 'spC', 'T', ...
           'diat_agg', 'diatFe', 'diatP', 'diatSi', 'diaz_agg', 'diaz_Nfix', 'diazFe', 'diazP', ...
            'photoFe_diat', 'photoFe_diaz', 'photoFe_sp', 'photoNO3_diat', 'photoNO3_diaz', ...
            'photoNO3_sp', 'PO4_diat_uptake', 'PO4_diaz_uptake', 'PO4_sp_uptake', 'sp_agg', ...
            'spFe', 'spP', 'zooC', 'photoC_NO3_TOT', 'photoC_TOT'}
        dim='4d'; %[x, y, z, t]
    case {'dustToSed', 'SSH', 'AODDUST', 'AODVIS', 'dst_a1SF', 'dst_a2SF', 'dst_a3SF', ...
            'PRECT', 'PSL', 'SST',  'BSF', 'HMXL', 'photoC_TOT_zint', 'photoC_TOT_zint_100m', ...
            'ATM_COARSE_DUST_FLUX_CPL', 'ATM_FINE_DUST_FLUX_CPL', 'SEAICE_DUST_FLUX_CPL', ...
           'P_iron_FLUX_100m', 'POC_FLUX_100m', 'POP_FLUX_100m', 'IRON_FLUX', ...
           'SiO2_FLUX_100m', 'TMXL_DR', 'XMXL_DR', ...
           'diatC_zint_100m_2', 'diazC_zint_100m_2', 'spC_zint_100m_2', ...
           'DSTFLXT', 'DSTDEP', 'DSL', 'QSOIL', 'TWS',...
           'TOTAL_DISCHARGE_TO_OCEAN_LIQ', ...
            'TOTAL_DISCHARGE_TO_OCEAN_ICE', 'NOx_FLUX', 'HBLT', 'HMXL_DR', 'TBLT', 'TMXL', 'XBLT', 'XMXL', ...
             'dry_deposition_NOy_as_N', 'dry_deposition_NHx_as_N', ...
            'dst_a1_SRF', 'dst_a2DDF', 'dst_a2SFWET', 'dst_a2_SRF',  ...
            'dst_a3_SRF', 'FSNS', 'FSDS', 'SFdst_a1', 'SFdst_a2', 'SFdst_a3', 'TAUX', 'TAUY', ...
            'TS', 'U10', 'wet_deposition_NHx_as_N', 'wet_deposition_NOy_as_N', ...
            'IVT', 'uIVT', 'vIVT', 'PRECL', 'PRECC', ...
            'PS', 'U850', 'UBOT', 'V850', 'VBOT', ...
            'diat_agg_zint_100m', 'diat_Fe_lim_Cweight_avg_100m', ...
            'diat_Fe_lim_surf', 'diat_light_lim_Cweight_avg_100m', 'diat_light_lim_surf', ...
            'diat_loss_zint_100m', 'diat_N_lim_Cweight_avg_100m', 'diat_N_lim_surf', ...
            'diat_P_lim_Cweight_avg_100m', 'diat_P_lim_surf', 'diat_SiO3_lim_Cweight_avg_100m', ...
            'diat_SiO3_lim_surf', 'diaz_agg_zint_100m', 'diaz_Fe_lim_Cweight_avg_100m', ...
            'diaz_Fe_lim_surf', 'diaz_light_lim_Cweight_avg_100m', 'diaz_light_lim_surf', ...
            'diaz_loss_zint', 'diaz_P_lim_Cweight_avg_100m', 'diaz_P_lim_surf', ...
            'graze_diat_zint_100m', 'graze_diat_zoo_zint_100m', 'graze_diaz_zint_100m', ...
            'graze_diaz_zoo_zint_100m', 'graze_sp_zint_100m', 'graze_sp_zoo_zint_100m', ...
            'LWUP_F', 'O2_ZMIN_DEPTH', 'photoC_diat_zint_100m', 'photoC_diaz_zint_100m', ...
            'photoC_NO3_TOT_zint_100m', 'photoC_sp_zint_100m', ...
            'sp_agg_zint_100m', 'sp_Fe_lim_Cweight_avg_100m', 'sp_Fe_lim_surf', ...
            'sp_light_lim_Cweight_avg_100m', 'sp_light_lim_surf', 'sp_loss_zint_100m', ...
            'sp_N_lim_Cweight_avg_100m', 'sp_N_lim_surf', 'sp_P_lim_Cweight_avg_100m', ...
            'sp_P_lim_surf', 'zoo_loss_zint_100m', 'NPP', 'GPP', 'TOTVEGC', 'aice', 'sithick'}
        dim='3d'; % [x, y, t]
end

end