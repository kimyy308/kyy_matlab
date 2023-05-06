clear all; close all; clc;
%% set path
[error_status, tmp.hostname] = system('hostname');
tmp.hostname=tmp.hostname(1:end-1);
switch tmp.hostname
    case 'Yong-Yubs-iMac-Pro.local'
        tmp.dropboxpath = '/Volumes/kyy_raid/kimyy/Dropbox';
    case {'da1', 'da2', 'da3', 'da4'}
        tmp.dropboxpath = '/mnt/lustre/proj/kimyy/Dropbox';
end
tmp.fs=filesep;
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);


%% year, month configuration
cfg.y=1970:2020;
cfg.m=1:12;
dirs.fdir_NINO34='/Volumes/kyy_raid/kimyy/Observation/EN4/NINO34';
dirs.fdir_NINO3='/Volumes/kyy_raid/kimyy/Observation/EN4/NINO3';
dirs.fdir_TROP='/Volumes/kyy_raid/kimyy/Observation/EN4/TROP';

a=0;
for ind_y=cfg.y
    str.y=num2str(ind_y,'%04i');
    for ind_m=cfg.m
        str.m=num2str(ind_m,'%02i');
        inds.t=(ind_y-min(cfg.y))*12+ind_m;
        name.file_NINO34_fldm=[dirs.fdir_NINO34, filesep, 'NINO34_EN4.2.2.fldmean.',str.y,str.m,'.nc'];
        data.nino34(inds.t)=ncread(name.file_NINO34_fldm,'temperature', [1 1 1 1], [inf, inf, 1, inf])-273.15;
        name.file_NINO34_temp=[dirs.fdir_NINO34, filesep, 'NINO34_EN4.2.2.',str.y,str.m,'.nc'];
        name.file_NINO3_temp=[dirs.fdir_NINO3, filesep, 'NINO3_EN4.2.2.',str.y,str.m,'.nc'];
        name.file_TROP_temp=[dirs.fdir_TROP, filesep, 'TROP_EN4.2.2.',str.y,str.m,'.nc'];
        
%         data.temp(:,:,:,inds.t)=ncread(name.file_NINO34_temp,'temperature', [1 1 1 1], [inf, inf, inf, inf])-273.15;
        data.temp(:,:,:,inds.t)=ncread(name.file_NINO3_temp,'temperature', [1 1 1 1], [inf, inf, inf, inf])-273.15;
        data.temp_TROP(:,:,:,inds.t)=ncread(name.file_TROP_temp,'temperature', [1 1 1 1], [inf, inf, inf, inf])-273.15;

        if ind_y==cfg.y(1) && ind_m==cfg.m(1)
            data.depth=ncread(name.file_NINO34_temp,'depth');
            data.depth_bnds=ncread(name.file_NINO34_temp,'depth_bnds');
            data.lon=ncread(name.file_NINO34_temp,'lon');
            data.lat=ncread(name.file_NINO34_temp,'lat');

            data.depth=ncread(name.file_NINO3_temp,'depth');
            data.depth_bnds=ncread(name.file_NINO3_temp,'depth_bnds');
            data.lon=ncread(name.file_NINO3_temp,'lon');
            data.lat=ncread(name.file_NINO3_temp,'lat');
            [data.lat_2d, data.lon_2d] = meshgrid(data.lat, data.lon);

            data.depth_TROP=ncread(name.file_TROP_temp,'depth');
            data.depth_bnds_TROP=ncread(name.file_TROP_temp,'depth_bnds');
            data.lon_TROP=ncread(name.file_TROP_temp,'lon');
            data.lat_TROP=ncread(name.file_TROP_temp,'lat');
            [data.lat_2d_TROP, data.lon_2d_TROP] = meshgrid(data.lat_TROP, data.lon_TROP);
        end
        for ind_lon=1:size(data.temp,1)
            for ind_lat=1:size(data.temp,2)
                tmp.t_1d=data.temp_TROP(ind_lon,ind_lat,:,inds.t);
                tmp.t_1d(tmp.t_1d<20)=NaN;
                if isnan(tmp.t_1d(1))
                    data.t_dep(ind_lon,ind_lat,inds.t)=NaN;
                    a=a+1;
%                     disp(['NaN ', num2str(a), 'ind_lon:', num2str(ind_lon), 'ind_lat:', num2str(ind_lat)]);
                else
                    data.t_dep(ind_lon,ind_lat,inds.t)=data.depth_bnds_TROP(2,find(isfinite(tmp.t_1d),1,'last'));
                end
            end
        end
        data.fldm_t_dep=Func_0011_get_area_weighted_mean(data.t_dep,data.lon_2d,data.lat_2d);

    end
end

data.nino34 = data.nino34 - mean(data.nino34); % make anomaly
data.AC_lag1=corr(data.nino34(1:end-1)', data.nino34(2:end)');
data.fldm_t_dep = data.fldm_t_dep - mean(data.fldm_t_dep); % make anomaly

data.dtdh=diff(data.fldm_t_dep); % get dtdh
data.dtdT=diff(data.nino34); % get dtdT
data.dtdh_resh=reshape(data.dtdh(1:end-11), [12, length(cfg.y)-1]); % reshape -> y,m
data.dtdT_resh=reshape(data.dtdT(1:end-11), [12, length(cfg.y)-1]); % reshape -> y,m

% plot(data.fldm_t_dep);
% corr(data.fldm_t_dep, data.nino34')o

corr(data.dtdh, data.nino34(1:end-1)')
scatter( data.nino34(1:end-1)', data.dtdh)
xlabel('T(^oC)'); ylabel('dH/dt(m/mon)');
% 
% 
% scatter( data.fldm_t_dep(1:end-1)', data.dtdh)
% xlabel('H(m)'); ylabel('dH/dt(m/mon)');

% scatter( data.nino34(1:end-1)', data.dtdT)
% xlabel('T(^oC)'); ylabel('dT/dt(^oC/mon)');
% 
% scatter( data.fldm_t_dep(1:end-1)', data.dtdT)
% xlabel('H(m)'); ylabel('dT/dt(^oC/mon)');

%% fitting R
% % % a=fit(data.fldm_t_dep(1:end-1), data.dtdT', 'poly1');
% % % coef.omega=a.p1;
a= fit(-data.nino34(1:end-1)', data.dtdh, 'poly1');
coef.R=a.p1;

%% fitting gamma as a function of month
data.nino34_resh=reshape(data.nino34(1:end), [12, length(cfg.y)]);
data.H_resh=reshape(data.fldm_t_dep(1:end), [12, length(cfg.y)]);
% for ind_m=cfg.m
%     data.dtdT_resh(ind_m, :) = diff(data.nino34_resh(ind_m,:));
%     data.dtdh_resh(ind_m, :) = diff(data.H_resh(ind_m,:));
%     data.var_season(ind_m)=var(data.nino34_resh(ind_m,:));
% end

% plot(data.var_season)
% xlabel('Mon'); ylabel('SST variance')

plot(data.dtdT_resh(1,:))
plot(data.dtdT_resh(6,:))

% % for ind_m=cfg.m
% %     a= fit(data.nino34_resh(ind_m,1:end-1)', data.dtdT_resh(ind_m,:)', 'poly1');
% %     coef.gamma_t(ind_m)=a.p1;
% % %     a= fit(data.nino34_resh(ind_m,1:end-1)', data.dtdh_resh(ind_m,:)', 'poly1');
% % %     coef.R_t(ind_m)=a.p1;
% %     a= fit(data.H_resh(ind_m,1:end-1)', data.dtdT_resh(ind_m,:)', 'poly1');
% %     coef.omega_t(ind_m)=a.p1;
% % end

%% multilinear regression (seasonal modulation)
for ind_m=cfg.m
    tmp.multi_set = [-2 .* data.nino34_resh(ind_m,1:end-1)' data.H_resh(ind_m,1:end-1)'];  %regress, -2T, H
    tmp.bb = regress(data.dtdT_resh(ind_m,:)', tmp.multi_set); % regress(Y, [x1 x2])
    coef.gamma_t(ind_m)=tmp.bb(1);
    coef.omega_t(ind_m)=tmp.bb(2);
end

%% multilinear regression (seasonal modulation x)
tmp.multi_set = [-2 .* data.nino34(1:end-1)', data.fldm_t_dep(1:end-1)];  %regress, -2T, H
tmp.bb = regress(data.dtdT(:), tmp.multi_set); % regress(Y, [x1 x2])
coef.gamma=tmp.bb(1);
coef.omega=tmp.bb(2);


plot(coef.gamma_t)
xlabel('Month'); ylabel('gamma(/yr)');
hold on
plot(1:12,repmat(mean(coef.gamma_t),[1 12]))
hold off

% plot(coef.R_t)
% xlabel('Month'); ylabel('R_t(/yr)');
% hold on
% plot(1:12,repmat(mean(coef.R_t),[1 12]))
% hold off

plot(coef.omega_t)
xlabel('Month'); ylabel('omega_t(/yr)');
hold on
plot(1:12,repmat(mean(coef.omega_t),[1 12]))
hold off



%% 3-year hindcast daily integration (seasonally modulated growth rate)
tmp.hc_y=3;
clear prog
for ind_y=cfg.y-min(cfg.y)+1
    tmp.y=ind_y+min(cfg.y)-1;
    

    %% Parametric recharge oscillator
    %% initialization
    prog.nino34(1:14)=NaN;
    prog.H(1:14)=NaN;
    prog.dT(1:14)=NaN;
    prog.dH(1:14)=NaN;
    prog.nino34(15)=data.nino34(((ind_y-1)*12)+1);
    prog.H(15)=data.fldm_t_dep(1);
    prog.noise=std(coef.gamma_t);
    prog.dt=1/30; % daily integration (360 days per year), but coefs are fitted to monthly data;
    
    %% seasonally modulated
    for ind_d=15:360*tmp.hc_y-1
        prog.ind_m(ind_d) = floor(mod(ind_d-0.1,360)/30)+1; % monthly index for seasonally modulating;
        prog.dT(ind_d) = -2* prog.dt.* (coef.gamma_t(prog.ind_m(ind_d)).*prog.nino34(ind_d)) ...
            + prog.dt.*coef.omega_t(prog.ind_m(ind_d)).*prog.H(ind_d) + prog.dt.*rand(1).*prog.noise;
        prog.dH(ind_d) = prog.dt.* (-coef.R.*prog.nino34(ind_d));
        prog.nino34(ind_d+1)= prog.nino34(ind_d) + prog.dT(ind_d);
        prog.H(ind_d+1)= prog.H(ind_d) + prog.dH(ind_d);
    end
    %% not seasonally modulated
    for ind_d=15:360*tmp.hc_y-1
        prog.ind_m(ind_d) = floor(mod(ind_d-0.1,360)/30)+1; % monthly index for seasonally modulating;
        prog.dT(ind_d) = -2* prog.dt.* (coef.gamma.*prog.nino34(ind_d)) ...
            + prog.dt.*coef.omega.*prog.H(ind_d) + prog.dt.*rand(1).*prog.noise;
        prog.dH(ind_d) = prog.dt.* (-coef.R.*prog.nino34(ind_d));
        prog.nino34(ind_d+1)= prog.nino34(ind_d) + prog.dT(ind_d);
        prog.H(ind_d+1)= prog.H(ind_d) + prog.dH(ind_d);
    end



    prog.t=prog.dt:prog.dt:360*length(cfg.y)*prog.dt;
    prog.ty=prog.t./12;
    prog.ty_mon=reshape(prog.ty, [30,length(prog.ty)/30]);
    prog.ty_mon=mean(prog.ty_mon,1,'omitnan');
    prog.nino34_mon=reshape(prog.nino34, [30,length(prog.nino34)/30]);
    prog.nino34_mon=mean(prog.nino34_mon,1,'omitnan');
    
    plot(prog.ty(1:360*tmp.hc_y), prog.nino34(1:360*tmp.hc_y))
    ylabel('Nino3.4(^oC)'); xlabel('Time(year)');
    grid minor
    name.modelvar=['hcst_nino34_i',num2str(tmp.y)];
    prog.(name.modelvar)=prog.nino34_mon;
    if ((ind_y-1)*12)+1+tmp.hc_y*12-1 > length(data.nino34)
        obs.(name.modelvar)(1:tmp.hc_y*12)=NaN;
    else
        obs.(name.modelvar)(1:tmp.hc_y*12)=data.nino34(((ind_y-1)*12)+1 : ((ind_y-1)*12)+tmp.hc_y*12);
    end
    %% assign variables according to lead month
    for lmonth=0:tmp.hc_y*12-1
        tmp.lmonth_str=num2str(lmonth, '%03i');
        name.modelvar_lm=['nino34m', '_model',  '_l', tmp.lmonth_str];
        tmp.mdata= prog.(name.modelvar)(lmonth+1);
        tmp.mind=(tmp.y-min(cfg.y))+1;
        prog.(name.modelvar_lm)(tmp.mind)=tmp.mdata;
%         prog.(name.modelvar_lm)(prog.(name.modelvar_lm)==0)=NaN;
        
        name.obsvar_lm=['nino34m', '_obs',  '_l', tmp.lmonth_str];
        tmp.mdata= obs.(name.modelvar)(lmonth+1);
        tmp.mind=(tmp.y-min(cfg.y))+1;
        obs.(name.obsvar_lm)(tmp.mind)=tmp.mdata;
%         data.(name.obsvar_lm)(data.(name.obsvar_lm)==0)=NaN;

        name.AR1var_lm=['nino34m', '_AR1',  '_l', tmp.lmonth_str];
        if lmonth==0
            prog.(name.AR1var_lm)(tmp.mind)=data.nino34(((ind_y-1)*12)+1);
        else
%             prog.(name.AR1var_lm)(tmp.mind)=Func_0030_AR1_prog(data.nino34(((ind_y-1)*12)+1), data.AC_lag1, lmonth, 0);
            prog.(name.AR1var_lm)(tmp.mind)=Func_0030_AR1_prog(data.nino34(((ind_y-1)*12)+1), data.AC_lag1, lmonth, prog.noise);
        end

    end
end
for lm=0:35
    tmp.lmonth_str=num2str(lm, '%03i');
    name.modelvar_lm=['nino34m', '_model',  '_l', tmp.lmonth_str];
    name.obsvar_lm=['nino34m', '_obs',  '_l', tmp.lmonth_str];
    tmp.ind=find(isfinite(obs.(name.obsvar_lm)));
    tmp.mod=prog.(name.modelvar_lm)(tmp.ind);
    tmp.obs=obs.(name.obsvar_lm)(tmp.ind);
    prog.corr_lm(lm+1)=corr(tmp.mod',tmp.obs');

    name.AR1var_lm=['nino34m', '_AR1',  '_l', tmp.lmonth_str];
    tmp.AR1=prog.(name.AR1var_lm)(tmp.ind);
    prog.corr_AR1_lm(lm+1)=corr(tmp.AR1',tmp.obs');
end
plot(0:35,prog.corr_lm(1:36), 'linewidth', 2)
hold on
plot(0:35,prog.corr_AR1_lm(1:36), 'linewidth', 2)
hold off
legend({'parametric recharge oscillator', 'AR1'})
ylabel('corr')
xlabel('lead month')
set(gca,'fontsize', 20)

plot(tmp.mod)
hold on
plot(tmp.AR1)
plot(tmp.obs)
hold off

mean((tmp.obs-mean(tmp.obs)).*(tmp.mod-mean(tmp.mod)))/(std(tmp.obs)*std(tmp.mod))
mean((tmp.AR1-mean(tmp.AR1)).*(tmp.obs-mean(tmp.obs)))/(std(tmp.AR1)*std(tmp.obs))


plot(prog.ty_mon(1:12*tmp.hc_y), prog.nino34_mon(1:12*tmp.hc_y))
ylabel('Nino3.4(^oC)'); xlabel('Time(year)');
grid minor
hold on
plot(prog.ty_mon(1:12*tmp.hc_y), data.nino34(1:12*tmp.hc_y))
hold off

% corr(prog.nino34_mon(1:12*5)', data.nino34(1:12*5)')

%% 3-year hindcast daily integration (fixed growth rate)
tmp.hc_y=3;
clear prog_fixG
for ind_y=cfg.y-min(cfg.y)+1
    tmp.y=ind_y+min(cfg.y)-1;
    %% initialization
    prog_fixG.nino34(1:14)=NaN;
    prog_fixG.H(1:14)=NaN;
    prog_fixG.dT(1:14)=NaN;
    prog_fixG.dH(1:14)=NaN;
    prog_fixG.nino34(15)=data.nino34(((ind_y-1)*12)+1);
    prog_fixG.H(15)=data.fldm_t_dep(1);
    
    prog_fixG.dt=1/30; % daily integration (360 days per year), but coefs are fitted to monthly data;
    
    for ind_d=15:360*tmp.hc_y-1
        prog_fixG.dT(ind_d) = -2* prog_fixG.dt.* (coef.gamma.*prog_fixG.nino34(ind_d) ...
            + prog_fixG.dt.*coef.omega.*prog_fixG.H(ind_d));
        prog_fixG.dH(ind_d) = prog_fixG.dt.* (-coef.R.*prog_fixG.nino34(ind_d));
        prog_fixG.nino34(ind_d+1)= prog_fixG.nino34(ind_d) + prog_fixG.dT(ind_d);
        prog_fixG.H(ind_d+1)= prog_fixG.H(ind_d) + prog_fixG.dH(ind_d);
    end
    prog_fixG.t=prog_fixG.dt:prog_fixG.dt:360*length(cfg.y)*prog_fixG.dt;
    prog_fixG.ty=prog_fixG.t./12;
    prog_fixG.ty_mon=reshape(prog_fixG.ty, [30,length(prog_fixG.ty)/30]);
    prog_fixG.ty_mon=mean(prog_fixG.ty_mon,1,'omitnan');
    prog_fixG.nino34_mon=reshape(prog_fixG.nino34, [30,length(prog_fixG.nino34)/30]);
    prog_fixG.nino34_mon=mean(prog_fixG.nino34_mon,1,'omitnan');
    
    plot(prog_fixG.ty(1:360*tmp.hc_y), prog_fixG.nino34(1:360*tmp.hc_y))
    ylabel('Nino3.4(^oC)'); xlabel('Time(year)');
    grid minor
    name.modelvar=['hcst_nino34_i',num2str(tmp.y)];
    prog_fixG.(name.modelvar)=prog_fixG.nino34_mon;
    if ((ind_y-1)*12)+1+tmp.hc_y*12-1 > length(data.nino34)
        obs.(name.modelvar)(1:tmp.hc_y*12)=NaN;
    else
        obs.(name.modelvar)(1:tmp.hc_y*12)=data.nino34(((ind_y-1)*12)+1 : ((ind_y-1)*12)+1+tmp.hc_y*12-1);
    end
    %% assign variables according to lead year
    for lmonth=0:tmp.hc_y*12-1
        tmp.lmonth_str=num2str(lmonth, '%03i');
        name.modelvar_lm=['nino34m', '_model',  '_l', tmp.lmonth_str];
        tmp.mdata= prog_fixG.(name.modelvar)(lmonth+1);
        tmp.mind=(tmp.y-min(cfg.y))+1;
        prog_fixG.(name.modelvar_lm)(tmp.mind)=tmp.mdata;
%         prog_fixG.(name.modelvar_lm)(prog_fixG.(name.modelvar_lm)==0)=NaN;
        
        name.obsvar_lm=['nino34m', '_obs',  '_l', tmp.lmonth_str];
        tmp.mdata= obs.(name.modelvar)(lmonth+1);
        tmp.mind=(tmp.y-min(cfg.y))+1;
        obs.(name.obsvar_lm)(tmp.mind)=tmp.mdata;
%         data.(name.obsvar_lm)(data.(name.obsvar_lm)==0)=NaN;
    end
end
for lm=0:35
    tmp.lmonth_str=num2str(lm, '%03i');
    name.modelvar_lm=['nino34m', '_model',  '_l', tmp.lmonth_str];
    name.obsvar_lm=['nino34m', '_obs',  '_l', tmp.lmonth_str];
    tmp.ind=find(isfinite(obs.(name.obsvar_lm)));
    tmp.mod=prog_fixG.(name.modelvar_lm)(tmp.ind);
    tmp.obs=obs.(name.obsvar_lm)(tmp.ind);
    prog_fixG.corr_lm(lm+1)=corr(tmp.mod',tmp.obs');
    tmp.mod2=prog.(name.modelvar_lm)(tmp.ind);
    tmp.corr_lm(lm+1)=corr(tmp.mod',tmp.mod2')


    plot(tmp.mod)
    hold on
    plot(tmp.mod2)
    hold off

end
plot(0:11,prog_fixG.corr_lm(1:12))
hold on
plot(0:11,prog.corr_lm(1:12))
hold off

plot(prog_fixG.ty_mon(1:12*tmp.hc_y), prog_fixG.nino34_mon(1:12*tmp.hc_y))
ylabel('Nino3.4(^oC)'); xlabel('Time(year)');
grid minor
hold on
plot(prog_fixG.ty_mon(1:12*tmp.hc_y), data.nino34(1:12*tmp.hc_y))
hold off








% % %% integration by obs value
% % for ind_y=cfg.y
% %     for ind_m=cfg.m
% %         inds.t=(ind_y-1950)*12+ind_m;
% %         prog_obs.dtdT(inds.t) = -2.*coef.gamma_t(ind_m).*data.nino34(inds.t) + coef.omega_t(ind_m).*data.fldm_t_dep(inds.t);
% %     end
% % end
% % corr(data.dtdT', prog_obs.dtdT(1:end-1)')

% plot(data.dtdT)
% hold on
% plot(prog_obs.dtdT(1:end-1))
% hold off

% data.temp size : 51 11 42 864

% tmptmp=data.temp(25,5,:,1);
% tmptmp(tmptmp<20)=NaN;
% find(isfinite(tmptmp), 1, 'last')
% squeeze(data.temp(25,5,:,1))