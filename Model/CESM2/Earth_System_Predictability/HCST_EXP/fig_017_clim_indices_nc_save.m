close all; clc; clear all;


% nccreate('CESM2_MP_sam.nc', ...
%                   'obs',...
%                   'Dimensions', {'time', inf}, ...
%                   'lens2',...
%                   'Dimensions', {'member_l', 50, 'time', inf}, ...;
%                   'assm',...
%                   'Dimensions', {'member_a', 20, 'time', inf}, ...
%                   'hcst',...
%                   'Dimensions', {'ly', 5, 'member_a', 20, 'time', inf}  );

ncdir='/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/nc';


date_time=[365/24: 365/12 :365/12*12*61];

name_var='SSH';
name_obs='CMEMS';

%% ENSO netcdf save
name_index='ENSO';

load(['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_', ...
    name_var, '_all_', name_index, '_obs_', name_obs,'.mat'])

fname_SST_ENSO=[ncdir, filesep, 'clim_indices_', name_var, '_all_', name_index, '_obs_', name_obs, '.nc'];
system(['rm ', fname_SST_ENSO]);
ncid = netcdf.create(fname_SST_ENSO,'NETCDF4');

%% define dimension
member_lens2_dimid = netcdf.defDim(ncid, 'member_LE', 50);
member_assm_dimid = netcdf.defDim(ncid,'member_ASSM', 20);
member_hcst_dimid = netcdf.defDim(ncid,'member_HCST', 20);
ly_dimid = netcdf.defDim(ncid,'ly', 5);
time_dimid = netcdf.defDim(ncid, 'time', 0);

%% global attribute
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                        'title', ['CESM2-MP ', 'climate indices']);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'used variable', cfg.var);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'domain', num2str(data_ENSO.regions));
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'used observation benchmark', cfg.obs_name);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'author', 'Created by Y.Y. Kim');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'date', date);

%% define variable
timevarid=netcdf.defVar(ncid, 'time', 'NC_DOUBLE', time_dimid);
    netcdf.putAtt(ncid,timevarid,'long_name','time');
    netcdf.putAtt(ncid,timevarid,'units','days since 1960-01-01 00:00:00');
    netcdf.putAtt(ncid,timevarid,'calendar','noleap');

member_lens2_varid=netcdf.defVar(ncid, 'members_LE', 'NC_STRING', member_lens2_dimid);
    netcdf.putAtt(ncid,member_lens2_varid,'long_name','CESM2 Large ensemble member list');
member_assm_varid=netcdf.defVar(ncid, 'members_ASSM', 'NC_STRING', member_assm_dimid);
    netcdf.putAtt(ncid,member_assm_varid,'long_name','CESM2 Assimilation(hindcast) member list');
member_hcst_varid=netcdf.defVar(ncid, 'members_HCST', 'NC_STRING', member_assm_dimid);
    netcdf.putAtt(ncid,member_hcst_varid,'long_name','CESM2 Hindcast member list');
ly_varid=netcdf.defVar(ncid, 'lead_year', 'NC_INT', ly_dimid);
    netcdf.putAtt(ncid,ly_varid,'long_name','Lead Year');

hcst_varid=netcdf.defVar(ncid, 'hcst', 'NC_FLOAT', [ly_dimid member_hcst_dimid time_dimid]);
                    netcdf.putAtt(ncid,hcst_varid,'long_name','hindcast');
assm_varid=netcdf.defVar(ncid, 'assm', 'NC_FLOAT', [member_assm_dimid time_dimid]);
                    netcdf.putAtt(ncid, assm_varid,'long_name','assimilation');
lens2_varid=netcdf.defVar(ncid, 'lens2', 'NC_FLOAT', [member_lens2_dimid time_dimid]);
                    netcdf.putAtt(ncid,lens2_varid,'long_name','large_ensemble');
obs_varid=netcdf.defVar(ncid, 'obs', 'NC_FLOAT', time_dimid);
                    netcdf.putAtt(ncid,obs_varid,'long_name','observation');

netcdf.endDef(ncid);

%% put variable
netcdf.putVar(ncid, member_lens2_varid, 0, cfg_lens2.len_mem, cfg_lens2.members);
netcdf.putVar(ncid, member_assm_varid, 0, cfg_assm.len_mem, cfg_assm.members);
netcdf.putVar(ncid, ly_varid, 0, 5, 1:5);
netcdf.putVar(ncid, timevarid, 0, length(date_time), date_time);
netcdf.putVar(ncid, obs_varid, 0, length(data_ENSO.obs), data_ENSO.obs);
netcdf.putVar(ncid, lens2_varid, [0 0], [size(data_ENSO.lens2)], data_ENSO.lens2);
netcdf.putVar(ncid, assm_varid, [0 0], [size(data_ENSO.assm)], data_ENSO.assm);
netcdf.putVar(ncid, hcst_varid, [0 0 0], [size(data_ENSO.hcst)], data_ENSO.hcst);

%% close file
netcdf.close(ncid);


%% AMO netcdf save
name_index='AMO';

load(['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_', ...
    name_var, '_all_', name_index, '_obs_', name_obs,'.mat'])

fname_SST_AMO=[ncdir, filesep, 'clim_indices_', name_var, '_all_', name_index, '_obs_', name_obs, '.nc'];
system(['rm ', fname_SST_AMO]);
ncid = netcdf.create(fname_SST_AMO,'NETCDF4');

%% define dimension
member_lens2_dimid = netcdf.defDim(ncid, 'member_LE', 50);
member_assm_dimid = netcdf.defDim(ncid,'member_ASSM', 20);
member_hcst_dimid = netcdf.defDim(ncid,'member_HCST', 20);
ly_dimid = netcdf.defDim(ncid,'ly', 5);
time_dimid = netcdf.defDim(ncid, 'time', 0);

%% global attribute
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                        'title', ['CESM2-MP ', 'climate indices']);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'used variable', cfg.var);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'domain_GLO', num2str(data_AMO.GLO_regions));
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'domain_ATL', num2str(data_AMO.ATL_regions));
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'used observation benchmark', cfg.obs_name);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'author', 'Created by Y.Y. Kim');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'date', date);

%% define variable
timevarid=netcdf.defVar(ncid, 'time', 'NC_DOUBLE', time_dimid);
    netcdf.putAtt(ncid,timevarid,'long_name','time');
    netcdf.putAtt(ncid,timevarid,'units','days since 1960-01-01 00:00:00');
    netcdf.putAtt(ncid,timevarid,'calendar','noleap');

member_lens2_varid=netcdf.defVar(ncid, 'members_LE', 'NC_STRING', member_lens2_dimid);
    netcdf.putAtt(ncid,member_lens2_varid,'long_name','CESM2 Large ensemble member list');
member_assm_varid=netcdf.defVar(ncid, 'members_ASSM', 'NC_STRING', member_assm_dimid);
    netcdf.putAtt(ncid,member_assm_varid,'long_name','CESM2 Assimilation(hindcast) member list');
member_hcst_varid=netcdf.defVar(ncid, 'members_HCST', 'NC_STRING', member_assm_dimid);
    netcdf.putAtt(ncid,member_hcst_varid,'long_name','CESM2 Hindcast member list');
ly_varid=netcdf.defVar(ncid, 'lead_year', 'NC_INT', ly_dimid);
    netcdf.putAtt(ncid,ly_varid,'long_name','Lead Year');

hcst_varid=netcdf.defVar(ncid, 'hcst', 'NC_FLOAT', [ly_dimid member_hcst_dimid time_dimid]);
                    netcdf.putAtt(ncid,hcst_varid,'long_name','hindcast');
assm_varid=netcdf.defVar(ncid, 'assm', 'NC_FLOAT', [member_assm_dimid time_dimid]);
                    netcdf.putAtt(ncid, assm_varid,'long_name','assimilation');
lens2_varid=netcdf.defVar(ncid, 'lens2', 'NC_FLOAT', [member_lens2_dimid time_dimid]);
                    netcdf.putAtt(ncid,lens2_varid,'long_name','large_ensemble');
obs_varid=netcdf.defVar(ncid, 'obs', 'NC_FLOAT', time_dimid);
                    netcdf.putAtt(ncid,obs_varid,'long_name','observation');

hcst_lp_varid=netcdf.defVar(ncid, 'hcst_lp', 'NC_FLOAT', [ly_dimid member_hcst_dimid time_dimid]);
                    netcdf.putAtt(ncid,hcst_lp_varid,'long_name','hindcast');
assm_lp_varid=netcdf.defVar(ncid, 'assm_lp', 'NC_FLOAT', [member_assm_dimid time_dimid]);
                    netcdf.putAtt(ncid, assm_lp_varid,'long_name','assimilation');
lens2_lp_varid=netcdf.defVar(ncid, 'lens2_lp', 'NC_FLOAT', [member_lens2_dimid time_dimid]);
                    netcdf.putAtt(ncid,lens2_lp_varid,'long_name','large_ensemble');
obs_lp_varid=netcdf.defVar(ncid, 'obs_lp', 'NC_FLOAT', time_dimid);
                    netcdf.putAtt(ncid,obs_lp_varid,'long_name','observation');

netcdf.endDef(ncid);

%% put variable
netcdf.putVar(ncid, member_lens2_varid, 0, cfg_lens2.len_mem, cfg_lens2.members);
netcdf.putVar(ncid, member_assm_varid, 0, cfg_assm.len_mem, cfg_assm.members);
netcdf.putVar(ncid, ly_varid, 0, 5, 1:5);
netcdf.putVar(ncid, timevarid, 0, length(date_time), date_time);
netcdf.putVar(ncid, obs_varid, 0, length(data_AMO.obs_dseason), data_AMO.obs_dseason);
netcdf.putVar(ncid, lens2_varid, [0 0], [size(data_AMO.lens2_dseason)], data_AMO.lens2_dseason);
netcdf.putVar(ncid, assm_varid, [0 0], [size(data_AMO.assm_dseason)], data_AMO.assm_dseason);
netcdf.putVar(ncid, hcst_varid, [0 0 0], [size(data_AMO.hcst_dseason)], data_AMO.hcst_dseason);
netcdf.putVar(ncid, obs_lp_varid, 0, length(data_AMO.obs_dseason_lp), data_AMO.obs_dseason_lp);
netcdf.putVar(ncid, lens2_lp_varid, [0 0], [size(data_AMO.lens2_dseason_lp)], data_AMO.lens2_dseason_lp);
netcdf.putVar(ncid, assm_lp_varid, [0 0], [size(data_AMO.assm_dseason_lp)], data_AMO.assm_dseason_lp);
netcdf.putVar(ncid, hcst_lp_varid, [0 0 0], [size(data_AMO.hcst_dseason_lp)], data_AMO.hcst_dseason_lp);

%% close file
netcdf.close(ncid);






%% PDO netcdf save
name_index='PDO';

load(['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_', ...
    name_var, '_all_', name_index, '_obs_', name_obs,'.mat'])

%% sign ordering
if strcmp(cfg.var,'SSH')
    nind=397;
else
    nind=1;
end

% ASSM
for mi=1:20
    tmp.corr=corrcoef(data_PDO.pct_obs(:,1), data_PDO.pct_assm(mi,nind:end,1), 'Rows', 'complete');
    if tmp.corr(1,2)<0
        data_PDO.pct_assm(mi,:,1)=-data_PDO.pct_assm(mi,:,1);
        data_PDO.lv_assm(mi,:,:,1)=-data_PDO.lv_assm(mi,:,:,1);
    end
end
% LENS2
for mi=1:50
    tmp.corr=corrcoef(data_PDO.pct_obs(:,1), data_PDO.pct_lens2(mi,nind:end,1), 'Rows', 'complete');
    if tmp.corr(1,2)<0
        data_PDO.pct_lens2(mi,:,1)=-data_PDO.pct_lens2(mi,:,1);
        data_PDO.lv_lens2(mi,:,:,1)=-data_PDO.lv_lens2(mi,:,:,1);
    end
end
% LENS2
for ly=1:5
for mi=1:20
    tmp.corr=corrcoef(data_PDO.pct_obs(:,1), data_PDO.pct_hcst(ly,mi,nind:end,1), 'Rows', 'complete');
    if tmp.corr(1,2)<0
        data_PDO.pct_hcst(ly,mi,:,1)=-data_PDO.pct_hcst(ly,mi,:,1);
        data_PDO.lv_hcst(ly,mi,:,:,1)=-data_PDO.lv_hcst(ly,mi,:,:,1);
    end
end
end

fname_SST_PDO=[ncdir, filesep, 'clim_indices_', name_var, '_all_', name_index, '_obs_', name_obs, '.nc'];
system(['rm ', fname_SST_PDO]);
ncid = netcdf.create(fname_SST_PDO,'NETCDF4');


%% define dimension
member_lens2_dimid = netcdf.defDim(ncid, 'member_LE', 50);
member_assm_dimid = netcdf.defDim(ncid,'member_ASSM', 20);
member_hcst_dimid = netcdf.defDim(ncid,'member_HCST', 20);
ly_dimid = netcdf.defDim(ncid,'ly', 5);
time_dimid = netcdf.defDim(ncid, 'time', 0);
one_dimid = netcdf.defDim(ncid, 'one', 1);
xi_dimid = netcdf.defDim(ncid, 'xi', data_PDO.cut_nlon);
yi_dimid = netcdf.defDim(ncid, 'yi', data_PDO.cut_nlat);

%% global attribute
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                        'title', ['CESM2-MP ', 'climate indices']);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'used variable', cfg.var);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'domain_GLO', num2str(data_AMO.GLO_regions));
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'domain_PAC', num2str(data_PDO.regions));
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'used observation benchmark', cfg.obs_name);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'author', 'Created by Y.Y. Kim');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'date', date);

%% define variable

xi_varid=netcdf.defVar(ncid, 'xi', 'NC_INT', xi_dimid);
                    netcdf.putAtt(ncid,xi_varid,'long_name','xi');
                    netcdf.putAtt(ncid,xi_varid,'axis','X');
yi_varid=netcdf.defVar(ncid, 'yi', 'NC_INT', yi_dimid);
                    netcdf.putAtt(ncid,yi_varid,'long_name','yi');
                    netcdf.putAtt(ncid,yi_varid,'axis','Y');

lon_varid=netcdf.defVar(ncid, 'Longitude', 'NC_FLOAT', [xi_dimid yi_dimid]);
    netcdf.putAtt(ncid,lon_varid,'long_name','Longitude');
lat_varid=netcdf.defVar(ncid, 'Latitude', 'NC_FLOAT', [xi_dimid yi_dimid]);
    netcdf.putAtt(ncid,lat_varid,'long_name','Latitude');

timevarid=netcdf.defVar(ncid, 'time', 'NC_DOUBLE', time_dimid);
    netcdf.putAtt(ncid,timevarid,'long_name','time');
    netcdf.putAtt(ncid,timevarid,'units','days since 1960-01-01 00:00:00');
    netcdf.putAtt(ncid,timevarid,'calendar','noleap');
    netcdf.putAtt(ncid,timevarid,'axis', 'T')

member_lens2_varid=netcdf.defVar(ncid, 'members_LE', 'NC_STRING', member_lens2_dimid);
    netcdf.putAtt(ncid,member_lens2_varid,'long_name','CESM2 Large ensemble member list');
    netcdf.putAtt(ncid,member_lens2_varid,'axis','M');
member_assm_varid=netcdf.defVar(ncid, 'members_ASSM', 'NC_STRING', member_assm_dimid);
    netcdf.putAtt(ncid,member_assm_varid,'long_name','CESM2 Assimilation(hindcast) member list');
    netcdf.putAtt(ncid,member_assm_varid,'axis','M');
member_hcst_varid=netcdf.defVar(ncid, 'members_HCST', 'NC_STRING', member_assm_dimid);
    netcdf.putAtt(ncid,member_hcst_varid,'long_name','CESM2 Hindcast member list');
    netcdf.putAtt(ncid,member_hcst_varid,'axis','M');
ly_varid=netcdf.defVar(ncid, 'lead_year', 'NC_INT', ly_dimid);
    netcdf.putAtt(ncid,ly_varid,'long_name','Lead Year');
    netcdf.putAtt(ncid,ly_varid,'axis','N');

hcst_varid=netcdf.defVar(ncid, 'hcst', 'NC_FLOAT', [time_dimid member_hcst_dimid ly_dimid]);
                    netcdf.putAtt(ncid,hcst_varid,'long_name','hindcast');
assm_varid=netcdf.defVar(ncid, 'assm', 'NC_FLOAT', [time_dimid member_assm_dimid ]);
                    netcdf.putAtt(ncid, assm_varid,'long_name','assimilation');
lens2_varid=netcdf.defVar(ncid, 'lens2', 'NC_FLOAT', [time_dimid member_lens2_dimid ]);
                    netcdf.putAtt(ncid,lens2_varid,'long_name','large_ensemble');
obs_varid=netcdf.defVar(ncid, 'obs', 'NC_FLOAT', time_dimid);
                    netcdf.putAtt(ncid,obs_varid,'long_name','observation');

hcst_lv_varid=netcdf.defVar(ncid, 'hcst_lv', 'NC_FLOAT', [xi_dimid yi_dimid member_hcst_dimid ly_dimid ]);
                    netcdf.putAtt(ncid,hcst_lv_varid,'long_name','hindcast');
assm_lv_varid=netcdf.defVar(ncid, 'assm_lv', 'NC_FLOAT', [xi_dimid yi_dimid member_assm_dimid ]);
                    netcdf.putAtt(ncid, assm_lv_varid,'long_name','assimilation');
lens2_lv_varid=netcdf.defVar(ncid, 'lens2_lv', 'NC_FLOAT', [xi_dimid yi_dimid member_lens2_dimid ]);
                    netcdf.putAtt(ncid,lens2_lv_varid,'long_name','large_ensemble');
obs_lv_varid=netcdf.defVar(ncid, 'obs_lv', 'NC_FLOAT', [xi_dimid yi_dimid]);
                    netcdf.putAtt(ncid,obs_lv_varid,'long_name','observation');




% hcst_lp_varid=netcdf.defVar(ncid, 'hcst_lp', 'NC_FLOAT', [ly_dimid member_hcst_dimid time_dimid]);
%                     netcdf.putAtt(ncid,hcst_lp_varid,'long_name','hindcast');
% assm_lp_varid=netcdf.defVar(ncid, 'assm_lp', 'NC_FLOAT', [member_assm_dimid time_dimid]);
%                     netcdf.putAtt(ncid, assm_lp_varid,'long_name','assimilation');
% lens2_lp_varid=netcdf.defVar(ncid, 'lens2_lp', 'NC_FLOAT', [member_lens2_dimid time_dimid]);
%                     netcdf.putAtt(ncid,lens2_lp_varid,'long_name','large_ensemble');
% obs_lp_varid=netcdf.defVar(ncid, 'obs_lp', 'NC_FLOAT', time_dimid);
%                     netcdf.putAtt(ncid,obs_lp_varid,'long_name','observation');

hcst_exp_varid=netcdf.defVar(ncid, 'hcst_exp', 'NC_FLOAT', [member_hcst_dimid ly_dimid ]);
                    netcdf.putAtt(ncid,hcst_exp_varid,'long_name','hindcast');
assm_exp_varid=netcdf.defVar(ncid, 'assm_exp', 'NC_FLOAT', member_assm_dimid);
                    netcdf.putAtt(ncid, assm_exp_varid,'long_name','assimilation');
lens2_exp_varid=netcdf.defVar(ncid, 'lens2_exp', 'NC_FLOAT', member_lens2_dimid);
                    netcdf.putAtt(ncid,lens2_exp_varid,'long_name','large_ensemble');
obs_exp_varid=netcdf.defVar(ncid, 'obs_exp', 'NC_FLOAT', one_dimid);
                    netcdf.putAtt(ncid,obs_exp_varid,'long_name','observation');

netcdf.endDef(ncid);

%% put variable
netcdf.putVar(ncid, member_lens2_varid, 0, cfg_lens2.len_mem, cfg_lens2.members);
netcdf.putVar(ncid, member_assm_varid, 0, cfg_assm.len_mem, cfg_assm.members);
netcdf.putVar(ncid, ly_varid, 0, 5, 1:5);
netcdf.putVar(ncid, timevarid, 0, length(date_time), date_time);
netcdf.putVar(ncid, lon_varid, [0 0], size(data_PDO.cut_tlong), data_PDO.cut_tlong);
netcdf.putVar(ncid, lat_varid, [0 0], size(data_PDO.cut_tlat), data_PDO.cut_tlat);

netcdf.putVar(ncid, obs_varid, 0, size(data_PDO.pct_obs,1), data_PDO.pct_obs(:,1));
% netcdf.putVar(ncid, lens2_varid, [0 0], [size(data_PDO.pct_lens2,1:2)], data_PDO.pct_lens2(:,:,1));
% netcdf.putVar(ncid, assm_varid, [0 0], [size(data_PDO.pct_assm,1:2)], data_PDO.pct_assm(:,:,1));
% netcdf.putVar(ncid, hcst_varid, [0 0 0], [size(data_PDO.pct_hcst,1:3)], data_PDO.pct_hcst(:,:,:,1));
netcdf.putVar(ncid, lens2_varid, [0 0], flip([size(data_PDO.pct_lens2,1:2)]), data_PDO.pct_lens2(:,:,1)');
netcdf.putVar(ncid, assm_varid, [0 0], flip([size(data_PDO.pct_assm,1:2)]), data_PDO.pct_assm(:,:,1)');
netcdf.putVar(ncid, hcst_varid, [0 0 0], [size(permute(data_PDO.pct_hcst, [3 2 1 4]),1:3)], permute(data_PDO.pct_hcst(:,:,:,1), [3 2 1]));

netcdf.putVar(ncid, obs_exp_varid, 0, 1, data_PDO.var_exp_obs(1));
netcdf.putVar(ncid, lens2_exp_varid, 0, [size(data_PDO.pct_lens2,1)], data_PDO.var_exp_lens2(:,1));
netcdf.putVar(ncid, assm_exp_varid, 0, [size(data_PDO.pct_assm,1)], data_PDO.var_exp_assm(:,1));
netcdf.putVar(ncid, hcst_exp_varid, [0 0], flip([size(data_PDO.pct_hcst,1:2)]), data_PDO.var_exp_hcst(:,:,1)');

netcdf.putVar(ncid, obs_lv_varid, [0 0], [size(data_PDO.lv_obs,1:2)], data_PDO.lv_obs(:,:,1));
netcdf.putVar(ncid, lens2_lv_varid, [0 0 0],[size(permute(data_PDO.lv_lens2, [2 3 1 4]),1:3)], permute(data_PDO.lv_lens2(:,:,:,1),[2 3 1]));
netcdf.putVar(ncid, assm_lv_varid, [0 0 0], [size(permute(data_PDO.lv_assm, [2 3 1 4]),1:3)], permute(data_PDO.lv_assm(:,:,:,1),[2 3 1]));
netcdf.putVar(ncid, hcst_lv_varid, [0 0 0 0], [size(permute(data_PDO.lv_hcst, [3 4 2 1 5]),1:4)], permute(data_PDO.lv_hcst(:,:,:,:,1),[3 4 2 1]));


%% close file
netcdf.close(ncid);


% abc=ncread(fname_SST_PDO, 'lens2');
% hold on
% for mi=1:50
%     plot(abc(:,mi))
% end







%% NPGO netcdf save
name_index='NPGO';

load(['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_', ...
    name_var, '_all_', 'PDO', '_obs_', name_obs,'.mat'])

nmode=2;
%% sign ordering
% ASSM
for mi=1:20
    tmp.corr=corrcoef(data_PDO.pct_obs(:,nmode), data_PDO.pct_assm(mi,nind:end,nmode), 'Rows', 'complete');
    if tmp.corr(1,2)<0
        data_PDO.pct_assm(mi,:,nmode)=-data_PDO.pct_assm(mi,:,nmode);
        data_PDO.lv_assm(mi,:,:,nmode)=-data_PDO.lv_assm(mi,:,:,nmode);
    end
end
% LENS2
for mi=1:50
    tmp.corr=corrcoef(data_PDO.pct_obs(:,nmode), data_PDO.pct_lens2(mi,nind:end,nmode), 'Rows', 'complete');
    if tmp.corr(1,2)<0
        data_PDO.pct_lens2(mi,:,nmode)=-data_PDO.pct_lens2(mi,:,nmode);
        data_PDO.lv_lens2(mi,:,:,nmode)=-data_PDO.lv_lens2(mi,:,:,nmode);
    end
end
% LENS2
for ly=1:5
for mi=1:20
    tmp.corr=corrcoef(data_PDO.pct_obs(:,nmode), data_PDO.pct_hcst(ly,mi,nind:end,nmode), 'Rows', 'complete');
    if tmp.corr(1,2)<0
        data_PDO.pct_hcst(ly,mi,:,nmode)=-data_PDO.pct_hcst(ly,mi,:,nmode);
        data_PDO.lv_hcst(ly,mi,:,:,nmode)=-data_PDO.lv_hcst(ly,mi,:,:,nmode);
    end
end
end

fname_SST_NPGO=[ncdir, filesep, 'clim_indices_', name_var, '_all_', name_index, '_obs_', name_obs, '.nc'];
system(['rm ', fname_SST_NPGO]);
ncid = netcdf.create(fname_SST_NPGO,'NETCDF4');


%% define dimension
member_lens2_dimid = netcdf.defDim(ncid, 'member_LE', 50);
member_assm_dimid = netcdf.defDim(ncid,'member_ASSM', 20);
member_hcst_dimid = netcdf.defDim(ncid,'member_HCST', 20);
ly_dimid = netcdf.defDim(ncid,'ly', 5);
time_dimid = netcdf.defDim(ncid, 'time', 0);
one_dimid = netcdf.defDim(ncid, 'one', 1);
xi_dimid = netcdf.defDim(ncid, 'xi', data_PDO.cut_nlon);
yi_dimid = netcdf.defDim(ncid, 'yi', data_PDO.cut_nlat);

%% global attribute
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                        'title', ['CESM2-MP ', 'climate indices']);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'used variable', cfg.var);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'domain_GLO', num2str(data_AMO.GLO_regions));
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'domain_PAC', num2str(data_PDO.regions));
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'used observation benchmark', cfg.obs_name);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'author', 'Created by Y.Y. Kim');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'date', date);

%% define variable

xi_varid=netcdf.defVar(ncid, 'xi', 'NC_INT', xi_dimid);
                    netcdf.putAtt(ncid,xi_varid,'long_name','xi');
                    netcdf.putAtt(ncid,xi_varid,'axis','X');
yi_varid=netcdf.defVar(ncid, 'yi', 'NC_INT', yi_dimid);
                    netcdf.putAtt(ncid,yi_varid,'long_name','yi');
                    netcdf.putAtt(ncid,yi_varid,'axis','Y');

lon_varid=netcdf.defVar(ncid, 'Longitude', 'NC_FLOAT', [xi_dimid yi_dimid]);
    netcdf.putAtt(ncid,lon_varid,'long_name','Longitude');
lat_varid=netcdf.defVar(ncid, 'Latitude', 'NC_FLOAT', [xi_dimid yi_dimid]);
    netcdf.putAtt(ncid,lat_varid,'long_name','Latitude');

timevarid=netcdf.defVar(ncid, 'time', 'NC_DOUBLE', time_dimid);
    netcdf.putAtt(ncid,timevarid,'long_name','time');
    netcdf.putAtt(ncid,timevarid,'units','days since 1960-01-01 00:00:00');
    netcdf.putAtt(ncid,timevarid,'calendar','noleap');
    netcdf.putAtt(ncid,timevarid,'axis', 'T')

member_lens2_varid=netcdf.defVar(ncid, 'members_LE', 'NC_STRING', member_lens2_dimid);
    netcdf.putAtt(ncid,member_lens2_varid,'long_name','CESM2 Large ensemble member list');
    netcdf.putAtt(ncid,member_lens2_varid,'axis','M');
member_assm_varid=netcdf.defVar(ncid, 'members_ASSM', 'NC_STRING', member_assm_dimid);
    netcdf.putAtt(ncid,member_assm_varid,'long_name','CESM2 Assimilation(hindcast) member list');
    netcdf.putAtt(ncid,member_assm_varid,'axis','M');
member_hcst_varid=netcdf.defVar(ncid, 'members_HCST', 'NC_STRING', member_assm_dimid);
    netcdf.putAtt(ncid,member_hcst_varid,'long_name','CESM2 Hindcast member list');
    netcdf.putAtt(ncid,member_hcst_varid,'axis','M');
ly_varid=netcdf.defVar(ncid, 'lead_year', 'NC_INT', ly_dimid);
    netcdf.putAtt(ncid,ly_varid,'long_name','Lead Year');
    netcdf.putAtt(ncid,ly_varid,'axis','N');

hcst_varid=netcdf.defVar(ncid, 'hcst', 'NC_FLOAT', [time_dimid member_hcst_dimid ly_dimid]);
                    netcdf.putAtt(ncid,hcst_varid,'long_name','hindcast');
assm_varid=netcdf.defVar(ncid, 'assm', 'NC_FLOAT', [time_dimid member_assm_dimid ]);
                    netcdf.putAtt(ncid, assm_varid,'long_name','assimilation');
lens2_varid=netcdf.defVar(ncid, 'lens2', 'NC_FLOAT', [time_dimid member_lens2_dimid ]);
                    netcdf.putAtt(ncid,lens2_varid,'long_name','large_ensemble');
obs_varid=netcdf.defVar(ncid, 'obs', 'NC_FLOAT', time_dimid);
                    netcdf.putAtt(ncid,obs_varid,'long_name','observation');

hcst_lv_varid=netcdf.defVar(ncid, 'hcst_lv', 'NC_FLOAT', [xi_dimid yi_dimid member_hcst_dimid ly_dimid ]);
                    netcdf.putAtt(ncid,hcst_lv_varid,'long_name','hindcast');
assm_lv_varid=netcdf.defVar(ncid, 'assm_lv', 'NC_FLOAT', [xi_dimid yi_dimid member_assm_dimid ]);
                    netcdf.putAtt(ncid, assm_lv_varid,'long_name','assimilation');
lens2_lv_varid=netcdf.defVar(ncid, 'lens2_lv', 'NC_FLOAT', [xi_dimid yi_dimid member_lens2_dimid ]);
                    netcdf.putAtt(ncid,lens2_lv_varid,'long_name','large_ensemble');
obs_lv_varid=netcdf.defVar(ncid, 'obs_lv', 'NC_FLOAT', [xi_dimid yi_dimid]);
                    netcdf.putAtt(ncid,obs_lv_varid,'long_name','observation');




% hcst_lp_varid=netcdf.defVar(ncid, 'hcst_lp', 'NC_FLOAT', [ly_dimid member_hcst_dimid time_dimid]);
%                     netcdf.putAtt(ncid,hcst_lp_varid,'long_name','hindcast');
% assm_lp_varid=netcdf.defVar(ncid, 'assm_lp', 'NC_FLOAT', [member_assm_dimid time_dimid]);
%                     netcdf.putAtt(ncid, assm_lp_varid,'long_name','assimilation');
% lens2_lp_varid=netcdf.defVar(ncid, 'lens2_lp', 'NC_FLOAT', [member_lens2_dimid time_dimid]);
%                     netcdf.putAtt(ncid,lens2_lp_varid,'long_name','large_ensemble');
% obs_lp_varid=netcdf.defVar(ncid, 'obs_lp', 'NC_FLOAT', time_dimid);
%                     netcdf.putAtt(ncid,obs_lp_varid,'long_name','observation');

hcst_exp_varid=netcdf.defVar(ncid, 'hcst_exp', 'NC_FLOAT', [member_hcst_dimid ly_dimid ]);
                    netcdf.putAtt(ncid,hcst_exp_varid,'long_name','hindcast');
assm_exp_varid=netcdf.defVar(ncid, 'assm_exp', 'NC_FLOAT', member_assm_dimid);
                    netcdf.putAtt(ncid, assm_exp_varid,'long_name','assimilation');
lens2_exp_varid=netcdf.defVar(ncid, 'lens2_exp', 'NC_FLOAT', member_lens2_dimid);
                    netcdf.putAtt(ncid,lens2_exp_varid,'long_name','large_ensemble');
obs_exp_varid=netcdf.defVar(ncid, 'obs_exp', 'NC_FLOAT', one_dimid);
                    netcdf.putAtt(ncid,obs_exp_varid,'long_name','observation');

netcdf.endDef(ncid);

%% put variable
netcdf.putVar(ncid, member_lens2_varid, 0, cfg_lens2.len_mem, cfg_lens2.members);
netcdf.putVar(ncid, member_assm_varid, 0, cfg_assm.len_mem, cfg_assm.members);
netcdf.putVar(ncid, ly_varid, 0, 5, 1:5);
netcdf.putVar(ncid, timevarid, 0, length(date_time), date_time);
netcdf.putVar(ncid, lon_varid, [0 0], size(data_PDO.cut_tlong), data_PDO.cut_tlong);
netcdf.putVar(ncid, lat_varid, [0 0], size(data_PDO.cut_tlat), data_PDO.cut_tlat);

netcdf.putVar(ncid, obs_varid, 0, size(data_PDO.pct_obs,1), data_PDO.pct_obs(:,nmode));
netcdf.putVar(ncid, lens2_varid, [0 0], flip([size(data_PDO.pct_lens2,1:2)]), data_PDO.pct_lens2(:,:,nmode)');
netcdf.putVar(ncid, assm_varid, [0 0], flip([size(data_PDO.pct_assm,1:2)]), data_PDO.pct_assm(:,:,nmode)');
netcdf.putVar(ncid, hcst_varid, [0 0 0], [size(permute(data_PDO.pct_hcst, [3 2 1 4]),1:3)], permute(data_PDO.pct_hcst(:,:,:,nmode), [3 2 1]));

netcdf.putVar(ncid, obs_exp_varid, 0, 1, data_PDO.var_exp_obs(nmode));
netcdf.putVar(ncid, lens2_exp_varid, 0, [size(data_PDO.pct_lens2,1)], data_PDO.var_exp_lens2(:,nmode));
netcdf.putVar(ncid, assm_exp_varid, 0, [size(data_PDO.pct_assm,1)], data_PDO.var_exp_assm(:,nmode));
netcdf.putVar(ncid, hcst_exp_varid, [0 0], flip([size(data_PDO.pct_hcst,1:2)]), data_PDO.var_exp_hcst(:,:,nmode)');

netcdf.putVar(ncid, obs_lv_varid, [0 0], [size(data_PDO.lv_obs,1:2)], data_PDO.lv_obs(:,:,nmode));
netcdf.putVar(ncid, lens2_lv_varid, [0 0 0],[size(permute(data_PDO.lv_lens2, [2 3 1 4]),1:3)], permute(data_PDO.lv_lens2(:,:,:,nmode),[2 3 1]));
netcdf.putVar(ncid, assm_lv_varid, [0 0 0], [size(permute(data_PDO.lv_assm, [2 3 1 4]),1:3)], permute(data_PDO.lv_assm(:,:,:,nmode),[2 3 1]));
netcdf.putVar(ncid, hcst_lv_varid, [0 0 0 0], [size(permute(data_PDO.lv_hcst, [3 4 2 1 5]),1:4)], permute(data_PDO.lv_hcst(:,:,:,:,nmode),[3 4 2 1]));


%% close file
netcdf.close(ncid);

figure;
abc=ncread(fname_SST_NPGO, 'assm');
hold on
for mi=1:20
    plot(data_PDO.time, abc(:,mi))
end


