close all; clear all; clc;

dir.lroot = ['~/Desktop/backup/Research/Postdoc/03_IBS/', ...
    '2022_predictability_assimilation_run/Lorenz_predictability_example/Model_src_kyy'];

dir.hcst = [dir.lroot, '/', 'ensPRED'];
dir.ctl = [dir.lroot, '/', 'CTL'];

cfg.ensnum=100;
cfg.t_ini_min=500;
% cfg.t_ini_max=95;
cfg.t_ini_max=2000; %999:default
cfg.pred_len=50;
cfg.t_ini_min_CTL=cfg.t_ini_min;
% cfg.t_ini_max_CTL=100;
cfg.t_ini_max_CTL=cfg.t_ini_max+cfg.pred_len;
cfg.long_range=500;

ext_flag=0;

fname=['/Users/kimyy/Desktop/backup/Research/Postdoc/03_IBS', ...
    '/2022_predictability_assimilation_run/Lorenz_predictability_example/', ...
    'lorenz_v2_temp.mat'];

if exist(fname)==0

%% read control run
for t_ini=cfg.t_ini_min_CTL:cfg.t_ini_max_CTL
    str.t_ini=num2str(t_ini, '%04i');
    dir.ctl_tmp=[dir.ctl, '/', 'rest'];
    fname_ctl=[dir.ctl_tmp, '/', 'restart.', str.t_ini, '.nc'];
%% raw
    ncid=netcdf.open(fname_ctl, 'NOWRITE');
    xvarid=netcdf.inqVarID(ncid,'X');
    yvarid=netcdf.inqVarID(ncid,'Y');
    zvarid=netcdf.inqVarID(ncid,'Z');
    CTL.X(t_ini-cfg.t_ini_min_CTL+1)=netcdf.getVar(ncid,xvarid);
    CTL.Y(t_ini-cfg.t_ini_min_CTL+1)=netcdf.getVar(ncid,yvarid);
    CTL.Z(t_ini-cfg.t_ini_min_CTL+1)=netcdf.getVar(ncid,zvarid);
%     CTL_X.(['p',str.t_ini])=netcdf.getVar(ncid,xvarid);
%     CTL_Y.(['p',str.t_ini])=netcdf.getVar(ncid,yvarid);
%     CTL_Z.(['p',str.t_ini])=netcdf.getVar(ncid,zvarid);
    disp(num2str(t_ini))
    netcdf.close(ncid);
end

%% read control for trajectory
for t_ini=cfg.t_ini_min_CTL:cfg.t_ini_max_CTL
    str.t_ini=num2str(t_ini, '%04i');
    dir.ctl_tmp=[dir.ctl, '/', 'out'];
    fname_ctl=[dir.ctl_tmp, '/', 'out.', str.t_ini, '.nc'];
%% raw
    ncid=netcdf.open(fname_ctl, 'NOWRITE');
    xvarid=netcdf.inqVarID(ncid,'X');
    yvarid=netcdf.inqVarID(ncid,'Y');
    zvarid=netcdf.inqVarID(ncid,'Z');
    period=(t_ini-cfg.t_ini_min_CTL)*10+1:(t_ini-cfg.t_ini_min_CTL+1)*10;
    CTL_traj.X(period)=netcdf.getVar(ncid,xvarid);
    CTL_traj.Y(period)=netcdf.getVar(ncid,yvarid);
    CTL_traj.Z(period)=netcdf.getVar(ncid,zvarid);
%     CTL_X.(['p',str.t_ini])=netcdf.getVar(ncid,xvarid);
%     CTL_Y.(['p',str.t_ini])=netcdf.getVar(ncid,yvarid);
%     CTL_Z.(['p',str.t_ini])=netcdf.getVar(ncid,zvarid);
    disp(num2str(t_ini))
    netcdf.close(ncid);
end

%% read ensemble predictions
for t_ini=cfg.t_ini_min:cfg.t_ini_max
    str.t_ini=num2str(t_ini, '%04i');
    for ensnum=1:cfg.ensnum
        str_ens=num2str(ensnum, '%02i');
%         disp('pred:', str.t_ini, 'ens:',str_ens);
        dir.hcst_tmp=[dir.hcst, '/', str.t_ini, '/', 'ens',str_ens, '/', 'rest'];
        for t_pred=1:cfg.pred_len
            str.t_hcst=num2str(t_ini+t_pred-1, '%04i');
            fname_hcst=[dir.hcst_tmp, '/', 'restart.', str.t_hcst, '.nc'];
            ncid=netcdf.open(fname_hcst, 'NOWRITE');
            xvarid=netcdf.inqVarID(ncid,'X');
            HCST.X.(['ens',str_ens]).(['i',str.t_ini])(t_pred)=netcdf.getVar(ncid,xvarid);
            netcdf.close(ncid);
            %% raw
%                 HCST.X.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=single(ncread(fname_hcst, 'X'));
        end
        disp([str.t_ini, ', ens : ', str_ens]);
    end
end

%% read first long prediction
t_ini=cfg.t_ini_min;
str.t_ini=num2str(t_ini, '%04i');
HCST_traj.X.(['i',str.t_ini])(1:cfg.ensnum, 1:(cfg.t_ini_max_CTL-cfg.t_ini_min_CTL)*10)=NaN;
HCST_traj.Y.(['i',str.t_ini])(1:cfg.ensnum, 1:(cfg.t_ini_max_CTL-cfg.t_ini_min_CTL)*10)=NaN;
HCST_traj.Z.(['i',str.t_ini])(1:cfg.ensnum, 1:(cfg.t_ini_max_CTL-cfg.t_ini_min_CTL)*10)=NaN;

for ensnum=1:cfg.ensnum
    str_ens=num2str(ensnum, '%02i');
    dir.hcst_tmp=[dir.hcst, '/', str.t_ini, '/', 'ens',str_ens, '/', 'out'];
    for t_pred=1:cfg.t_ini_max_CTL-cfg.t_ini_min_CTL
        str.t_hcst=num2str(t_ini+t_pred-1, '%04i');
        fname_hcst=[dir.hcst_tmp, '/', 'out.', str.t_hcst, '.nc'];
        ncid=netcdf.open(fname_hcst, 'NOWRITE');
        xvarid=netcdf.inqVarID(ncid,'X');
        HCST_traj.X.(['i',str.t_ini])(ensnum,(t_pred-1)*10+1:(t_pred)*10)=netcdf.getVar(ncid,xvarid);
        yvarid=netcdf.inqVarID(ncid,'Y');
        HCST_traj.Y.(['i',str.t_ini])(ensnum,(t_pred-1)*10+1:(t_pred)*10)=netcdf.getVar(ncid,yvarid);
        zvarid=netcdf.inqVarID(ncid,'Z');
        HCST_traj.Z.(['i',str.t_ini])(ensnum,(t_pred-1)*10+1:(t_pred)*10)=netcdf.getVar(ncid,zvarid);
        netcdf.close(ncid);
    end
    disp([str.t_ini, ', ens : ', str_ens]);
end
HCST_traj.X_em.(['i',str.t_ini])=mean(HCST_traj.X.(['i',str.t_ini]),1);
HCST_traj.Y_em.(['i',str.t_ini])=mean(HCST_traj.Y.(['i',str.t_ini]),1);
HCST_traj.Z_em.(['i',str.t_ini])=mean(HCST_traj.Z.(['i',str.t_ini]),1);


    save(fname, '-v7.3');
else
    load(fname);
end







% %% data combine for CTL
% for t_ini=cfg.t_ini_min_CTL:cfg.t_ini_max_CTL
%     str.t_ini=num2str(t_ini, '%04i');
%     t_plot=t_ini:0.01:t_ini+0.99;
%     plot(t_plot, CTL_Z.(['p',str.t_ini]), 'b');
% 
%     data_comb.t((t_ini-cfg.t_ini_min_CTL)*100+1:(t_ini-cfg.t_ini_min_CTL+1)*100)=t_plot;
%     CTL.X((t_ini-cfg.t_ini_min_CTL)*100+1:(t_ini-cfg.t_ini_min_CTL+1)*100)=CTL_X.(['p',str.t_ini]);
%     data_comb.Y((t_ini-cfg.t_ini_min_CTL)*100+1:(t_ini-cfg.t_ini_min_CTL+1)*100)=CTL_Y.(['p',str.t_ini]);
%     data_comb.Z((t_ini-cfg.t_ini_min_CTL)*100+1:(t_ini-cfg.t_ini_min_CTL+1)*100)=CTL_Z.(['p',str.t_ini]);
% end
% hold off
% 
% for t_ini=cfg.t_ini_min:cfg.pred_len:cfg.t_ini_max
%     str.t_ini=num2str(t_ini, '%04i');
%     for ensi=1:cfg.ensnum
%         str.ens=num2str(ensi, '%02i');
%         for ly=1:cfg.pred_len
%             str.t_pred=num2str(t_ini+ly-1,'%04i');
%             hcst_comb.X.(['ens',str.ens]).(['i',str.t_ini])((ly-1)*100+1:(ly-1)*100+100) = ...
%                 HCST.X.(['ens',str.ens]).(['i',str.t_ini]).(['p',str.t_pred]);
%             hcst_comb_unini.X.(['ens',str.ens]).(['i',str.t_ini])((ly-1)*100+1:(ly-1)*100+100) = ...
%                 HCST_unini_X.(['ens',str.ens]).(['i',str.t_ini]).(['p',str.t_pred]);
%         end
%     end
% end

%% assign data by lead time
for ensi=1:cfg.ensnum
    str.ens=num2str(ensi, '%02i');
    for lt_i=1:50
        for t_ini=cfg.t_ini_min:cfg.t_ini_max
            str.t_ini=num2str(t_ini, '%04i');
            HCST.X.data_lt(ensi,lt_i,t_ini-cfg.t_ini_min+1) = ...
                HCST.X.(['ens',str.ens]).(['i',str.t_ini])(lt_i);
        end
    end
end
HCST.X.data_lt_em=squeeze(mean(HCST.X.data_lt,1));


%% corr coeficient calculation -------------------------------------------------------

% corrcoef(CTL.X(901:100:3401),  HCST.X.data_lt(1,1,6:31))  % (lt=1, 60~85)
% corrcoef(CTL.X(900:100:3400),  HCST.X.data_lt_em(500,1:26))  % (lt=500, 60~85)

% for HCST
tis=0:0.1:4.9;
for tti=1:length(tis)
    ti=tis(tti);
%     tmp.tmin=tti;
%     tmp.tmax=tti+499;
    tmp.tmin=tti;
    tmp.tmax=tti+(cfg.t_ini_max-cfg.t_ini_min);
    tmp.ts=1; % 1 + 10 ~ 1 + 11 - 10
    tmp.te=cfg.t_ini_max-cfg.t_ini_min+1; % 95 - 55 + 1 ~ 95 -55 + 1 - 10

    [tmp.corr, tmp.corr_p]=corrcoef(CTL.X(tmp.tmin:tmp.tmax),  ...
        HCST.X.data_lt_em(tti,tmp.ts:tmp.te)); 
    HCST.X.corr_em(tti)=tmp.corr(1,2);
    HCST.X.corr_em_p(tti)=tmp.corr_p(1,2);
    HCST.X.sig_t_em(tti)=std(HCST.X.data_lt_em(tti,tmp.ts:tmp.te));
    tmp.cov=cov(CTL.X(tmp.tmin:tmp.tmax),  ...
        squeeze(HCST.X.data_lt_em(tti,tmp.ts:tmp.te)));  
    HCST.X.cov_em(tti)=tmp.cov(1,2);
    for ensi=1:cfg.ensnum
        [tmp.corr, tmp.corr_p]=corrcoef(CTL.X(tmp.tmin:tmp.tmax),  ...
            HCST.X.data_lt(ensi,tti,tmp.ts:tmp.te));  
        HCST.X.corr(ensi,tti)=tmp.corr(1,2);
        HCST.X.corr_p(ensi,tti)=tmp.corr_p(1,2);
        HCST.X.sig_t_indv(ensi,tti)=std(HCST.X.data_lt(ensi,tti,tmp.ts:tmp.te));
        tmp.cov=cov(CTL.X(tmp.tmin:tmp.tmax),  ...
            squeeze(HCST.X.data_lt(ensi,tti,tmp.ts:tmp.te)));  
        HCST.X.cov_indv(ensi,tti)=tmp.cov(1,2);
    end
    HCST.X.corr_med(tti)=median(HCST.X.corr(:,tti));
    HCST.X.corr_mean(tti)=mean(HCST.X.corr(:,tti));
    HCST.X.sig_t_med(tti)=squeeze(median(HCST.X.sig_t_indv(:,tti),1));
    HCST.X.sig_t_mean(tti)=squeeze(mean(HCST.X.sig_t_indv(:,tti),1));
    HCST.X.cov_med(tti)=median(HCST.X.cov_indv(:,tti));
    HCST.X.cov_mean(tti)=mean(HCST.X.cov_indv(:,tti),1);
end


%% Brier score calculation
tis=0:0.1:4.9;
clear bs
for tti=1:length(tis) % lead time
    ti=tis(tti);
    tmp.tmin=tti;
    tmp.tmax=tti+(cfg.t_ini_max-cfg.t_ini_min);
    tmp.ts=1; % 1 + 10 ~ 1 + 11 - 10
    tmp.te=cfg.t_ini_max-cfg.t_ini_min+1; % 95 - 55 + 1 ~ 95 -55 + 1 - 10

    tmp.corr=corrcoef(CTL.X(tmp.tmin:tmp.tmax), HCST.X.data_lt(1,tti,:));
    bsf(tti)=tmp.corr(1,2);
    for fi=tmp.tmin:tmp.tmax % initialized time
%         sign(CTL.X(fi));
        tmp.fv=sign(HCST.X.data_lt(:,tti,fi-tmp.tmin+1));
%         tmp.fv=sign(HCST.X.data_lt(1,tti,fi-tti+1));
%% true rate
        tmp.fv2=tmp.fv;
        tmp.fv2(tmp.fv==sign(CTL.X(fi)))=1;
        tmp.fv2(tmp.fv~=sign(CTL.X(fi)))=0;
%         bs(tti,fi-tti+1)=tmp.fv2;
        
%% false rate
        tmp.fv3=tmp.fv;
        tmp.fv3(tmp.fv==sign(CTL.X(fi)))=0;
        tmp.fv3(tmp.fv~=sign(CTL.X(fi)))=1;

        bs(tti,fi-tti+1)=(mean(tmp.fv2)-1).^2 + (mean(tmp.fv3)-0).^2;
    end
    bsm(tti)=mean(bs(tti,:));
% end
% hold off
% plot(mean(bs,2));
% hold on
% plot(HCST.X.corr_em);
% plot(HCST.X.corr_med);

    [tmp.corr, tmp.corr_p]=corrcoef(CTL.X(tmp.tmin:tmp.tmax),  ...
        HCST.X.data_lt_em(tti,tmp.ts:tmp.te)); 
    HCST.X.corr_em(tti)=tmp.corr(1,2);
    HCST.X.corr_em_p(tti)=tmp.corr_p(1,2);
    HCST.X.sig_t_em(tti)=std(HCST.X.data_lt_em(tti,tmp.ts:tmp.te));
    tmp.cov=cov(CTL.X(tmp.tmin:tmp.tmax),  ...
        squeeze(HCST.X.data_lt_em(tti,tmp.ts:tmp.te)));  
    HCST.X.cov_em(tti)=tmp.cov(1,2);
    for ensi=1:cfg.ensnum
        [tmp.corr, tmp.corr_p]=corrcoef(CTL.X(tmp.tmin:tmp.tmax),  ...
            HCST.X.data_lt(ensi,tti,tmp.ts:tmp.te));  
        HCST.X.corr(ensi,tti)=tmp.corr(1,2);
        HCST.X.corr_p(ensi,tti)=tmp.corr_p(1,2);
        HCST.X.sig_t_indv(ensi,tti)=std(HCST.X.data_lt(ensi,tti,tmp.ts:tmp.te));
        tmp.cov=cov(CTL.X(tmp.tmin:tmp.tmax),  ...
            squeeze(HCST.X.data_lt(ensi,tti,tmp.ts:tmp.te)));  
        HCST.X.cov_indv(ensi,tti)=tmp.cov(1,2);
    end
    HCST.X.corr_med(tti)=median(HCST.X.corr(:,tti));
    HCST.X.corr_mean(tti)=mean(HCST.X.corr(:,tti));

end


%% corr coeficient calculation as a function of ensemble members -----------

% for HCST
% tis=0:0.1:4.9;
tis=2.0;
% for tti=1:length(tis)
    tti=21;
%     ti=tis(tti);
    ti=tis;
%     tmp.tmin=tti;
%     tmp.tmax=tti+499;
    tmp.tmin=tti;
    tmp.tmax=tti+(cfg.t_ini_max-cfg.t_ini_min);
    tmp.ts=1; % 1 + 10 ~ 1 + 11 - 10
    tmp.te=cfg.t_ini_max-cfg.t_ini_min+1; % 95 - 55 + 1 ~ 95 -55 + 1 - 10
    
    for mi=1:100
        for triali=1:100
            samples=datasample(1:100, mi, 'Replace',false);
            tmp.lt_em=mean(HCST.X.data_lt(samples,tti,tmp.ts:tmp.te),1);
            [tmp.corr, tmp.corr_p]=corrcoef(CTL.X(tmp.tmin:tmp.tmax),  ...
                tmp.lt_em);
            skills_em(mi,triali)=tmp.corr(1,2);
            
            skills_ri=NaN(1,mi);
            for mmi=1:mi
                [tmp.corr, tmp.corr_p]=corrcoef(CTL.X(tmp.tmin:tmp.tmax),  ...
                HCST.X.data_lt(samples(mmi),tti,tmp.ts:tmp.te));
                skills_ri(mmi)=tmp.corr(1,2);
            end
            skills_ri_med(mi,triali)=median(skills_ri);
            skills_ri_mean(mi,triali)=mean(skills_ri);
        end
        mi
    end
    
%%     box plots
    hold on
    for mi=1:100
        for triali=1:100
            scatter(mi,skills_ri_med(mi,triali));
        end
    end
    boxplot(skills_em', 1:100, 'BoxStyle', 'outline', 'symbol', '', 'Whisker', 2)
    hold on
    boxplot(skills_ri_med', 1:100, 'BoxStyle', 'outline', 'symbol', '', 'Whisker', 2)
    xlabel('number of ensemble members');
    ylabel('skill')



% % %     [tmp.corr, tmp.corr_p]=corrcoef(CTL.X(tmp.tmin:tmp.tmax),  ...
% % %         HCST.X.data_lt_em(tti,tmp.ts:tmp.te)); 
% % %     HCST.X.corr_em(tti)=tmp.corr(1,2);
% % %     HCST.X.corr_em_p(tti)=tmp.corr_p(1,2);
% % %     HCST.X.sig_t_em(tti)=std(HCST.X.data_lt_em(tti,tmp.ts:tmp.te));
% % %     tmp.cov=cov(CTL.X(tmp.tmin:tmp.tmax),  ...
% % %         squeeze(HCST.X.data_lt_em(tti,tmp.ts:tmp.te)));  
% % %     HCST.X.cov_em(tti)=tmp.cov(1,2);
% % %     for ensi=1:cfg.ensnum
% % %         [tmp.corr, tmp.corr_p]=corrcoef(CTL.X(tmp.tmin:tmp.tmax),  ...
% % %             HCST.X.data_lt(ensi,tti,tmp.ts:tmp.te));  
% % %         HCST.X.corr(ensi,tti)=tmp.corr(1,2);
% % %         HCST.X.corr_p(ensi,tti)=tmp.corr_p(1,2);
% % %         HCST.X.sig_t_indv(ensi,tti)=std(HCST.X.data_lt(ensi,tti,tmp.ts:tmp.te));
% % %         tmp.cov=cov(CTL.X(tmp.tmin:tmp.tmax),  ...
% % %             squeeze(HCST.X.data_lt(ensi,tti,tmp.ts:tmp.te)));  
% % %         HCST.X.cov_indv(ensi,tti)=tmp.cov(1,2);
% % %     end
% % %     HCST.X.corr_med(tti)=median(HCST.X.corr(:,tti));
% % %     HCST.X.corr_mean(tti)=mean(HCST.X.corr(:,tti));
% % %     HCST.X.sig_t_med(tti)=squeeze(median(HCST.X.sig_t_indv(:,tti),1));
% % %     HCST.X.sig_t_mean(tti)=squeeze(mean(HCST.X.sig_t_indv(:,tti),1));
% % %     HCST.X.cov_med(tti)=median(HCST.X.cov_indv(:,tti));
% % %     HCST.X.cov_mean(tti)=mean(HCST.X.cov_indv(:,tti),1);
% % % % end







%% figure start


fig_cfg.fig_size = [0,0,19,14.5]; %% paper size (original)
fig_cfg.cb_size = [1, 13, 17, 0.3];
loc_column_first=1;
loc_row_first=7;
fig_h = figure('name', 'fig_lorenz','PaperUnits','inches', ...
        'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','on');

%% plot using surf (CTL, OBS)
ax_m_1=subplot(7,7,1);
fig_cfg.ax_size = [loc_column_first, loc_row_first, 5, 5];
period=1:5000;
fig_weig=1;
% ax_m_1=axes('fontsize',14, 'fontname','freeserif'); 
set(ax_m_1, 'Parent', fig_h)
set(ax_m_1,'Units','inches','Position',fig_cfg.ax_size);
set(ax_m_1,'fontsize',12);
grid on
hold on
plots.X=CTL_traj.X(period);
plots.Y=CTL_traj.Y(period);
plots.Z=CTL_traj.Z(period);
tmp.parula=parula(200+(max(period))); % from light blue, until darker yellow
if isfield(plots, 'C0')
    plots=rmfield(plots,'C0');
end
for ti=1:length(plots.X)
    plots.C0(ti,1:2,1)=tmp.parula(ti+100,1);
    plots.C0(ti,1:2,2)=tmp.parula(ti+100,2);
    plots.C0(ti,1:2,3)=tmp.parula(ti+100,3);
end
 objsurf=surf(ax_m_1, [plots.X(:) plots.X(:)], [plots.Y(:) plots.Y(:)], [plots.Z(:) plots.Z(:)], ...
        plots.C0, ...  % Reshape and replicate data
     'FaceColor', 'none', ...    % Don't bother filling faces with color
     'EdgeColor', 'interp', ...  % Use interpolated color for edges
     'LineWidth', fig_weig, 'EdgeAlpha', 1);            % Make a thicker line
title1=title('(a) Observation', 'fontsize', 20);
view(210,23)
caxis([0 max(period)/100]);
xlabel('$$ x(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
xlim([-27 27])
ylabel('$$ y(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
ylim([-29 29]);
zlabel('$$ z(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
zlim([0 55]);
cb = colorbar(ax_m_1,'units','inches','Location', 'northoutside', 'position',fig_cfg.cb_size);
set(cb,'fontsize',15,'fontname','freeserif','TickDir','both');
cb_title=title(cb,'$$ \tau $$','fontsize', 22, 'Position', [1240, 0, 0]); % hor, ver, ?
set(cb_title, 'interpreter', 'latex');

%scatter
ax_m_2=subplot(7,7,2);
set(ax_m_2, 'Color','none');
set(ax_m_2, 'Parent', fig_h)
set(ax_m_2,'Units','inches','Position',fig_cfg.ax_size);
set(ax_m_2,'fontsize',12);
hold on
sc1=scatter3(ax_m_2,plots.X(1), plots.Y(1), plots.Z(1), 'MarkerFaceColor', plots.C0(1,1,:), ...
    'MarkerEdgeColor', 'none', 'LineWidth', 3);
view(210,23)
caxis([0 max(period)/100]);
xlabel('$$ x(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
xlim([-27 27])
ylabel('$$ y(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
ylim([-29 29]);
zlabel('$$ z(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
zlim([0 55]);
cb = colorbar(ax_m_2,'units','inches','Location', 'northoutside', 'position',fig_cfg.cb_size);
set(cb,'fontsize',15,'fontname','freeserif','TickDir','both');
cb_title=title(cb,'$$ \tau $$','fontsize', 22, 'Position', [1240, 0, 0]); % hor, ver, ?
set(cb_title, 'interpreter', 'latex');





%% plot using surf, gray noise (ensemble mean)
ax_m_3=subplot(7,7,3);
fig_cfg.ax_size = [loc_column_first + 6, loc_row_first, 5, 5];
period=1:5000;
fig_weig=1;
% ax_m_3=axes('fontsize',14, 'fontname','freeserif'); 
set(ax_m_3, 'Parent', fig_h)
set(ax_m_3,'Units','inches','Position',fig_cfg.ax_size);
set(ax_m_3,'fontsize',12);
hold on
grid on
hold all

for ensi=1:cfg.ensnum
% for ensi=1:30

    strens=num2str(ensi, '%02i');
    plots.X=HCST_traj.X.i0500(ensi,period);
    plots.Y=HCST_traj.Y.i0500(ensi,period);
    plots.Z=HCST_traj.Z.i0500(ensi,period);
    for ti=1:length(plots.X)
        plots.C0(ti,1:2,1)=0.9;
        plots.C0(ti,1:2,2)=0.9;
        plots.C0(ti,1:2,3)=0.9;
    end
    hold on;
    objsurf=surf(ax_m_3, [plots.X(:) plots.X(:)], [plots.Y(:) plots.Y(:)], [plots.Z(:) plots.Z(:)], ...
        plots.C0, ...  % Reshape and replicate data
     'FaceColor', 'none', ...    % Don't bother filling faces with color
     'EdgeColor', 'interp', ...  % Use interpolated color for edges
     'LineWidth', fig_weig, 'EdgeAlpha', 0.04);            % Make a thicker line
%     drawnow;
end
view(210,23)
title1=title('(b) Ensemble Mean', 'fontsize', 20);
sc1=scatter3(ax_m_3,plots.X(1), plots.Y(1), plots.Z(1), 'MarkerFaceColor', plots.C0(1,1,:), ...
    'MarkerEdgeColor', 'none', 'LineWidth', 3);
view(210,23)
caxis([0 max(period)/100]);
xlabel('$$ x(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
xlim([-27 27])
ylabel('$$ y(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
ylim([-29 29]);
zlabel('$$ z(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
zlim([0 55]);
cb = colorbar(ax_m_3,'units','inches','Location', 'northoutside', 'position',fig_cfg.cb_size);
set(cb,'fontsize',15,'fontname','freeserif','TickDir','both');
cb_title=title(cb,'$$ \tau $$','fontsize', 22, 'Position', [1240, 0, 0]); % hor, ver, ?
set(cb_title, 'interpreter', 'latex');


%% color ensemble mean
ax_m_4=subplot(7,7,4);
set(ax_m_4, 'Color','none');
fig_cfg.ax_size = [loc_column_first + 6, loc_row_first, 5, 5];
period=1:5000;
fig_weig=1;
% ax_m_4=axes('fontsize',14, 'fontname','freeserif'); 
set(ax_m_4, 'Parent', fig_h)
set(ax_m_4,'Units','inches','Position',fig_cfg.ax_size);
set(ax_m_4,'fontsize',12);
hold on
plots.X=HCST_traj.X_em.i0500(period);
plots.Y=HCST_traj.Y_em.i0500(period);
plots.Z=HCST_traj.Z_em.i0500(period);
tmp.parula=parula(200+(max(period))); % from light blue, until darker yellow
plots=rmfield(plots,'C0');
for ti=1:length(plots.X)
    plots.C0(ti,1:2,1)=tmp.parula(ti+100,1);
    plots.C0(ti,1:2,2)=tmp.parula(ti+100,2);
    plots.C0(ti,1:2,3)=tmp.parula(ti+100,3);
end
 objsurf=surf(ax_m_4, [plots.X(:) plots.X(:)], [plots.Y(:) plots.Y(:)], [plots.Z(:) plots.Z(:)], ...
        plots.C0, ...  % Reshape and replicate data
     'FaceColor', 'none', ...    % Don't bother filling faces with color
     'EdgeColor', 'interp', ...  % Use interpolated color for edges
     'LineWidth', fig_weig, 'EdgeAlpha', 1);            % Make a thicker line


title1=title('(b) Ensemble Mean', 'fontsize', 20);
sc1=scatter3(ax_m_4,plots.X(1), plots.Y(1), plots.Z(1), 'MarkerFaceColor', plots.C0(1,1,:), ...
    'MarkerEdgeColor', 'none', 'LineWidth', 3);
view(210,23)
caxis([0 max(period)/100]);
xlabel('$$ x(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
xlim([-27 27])
ylabel('$$ y(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
ylim([-29 29]);
zlabel('$$ z(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
zlim([0 55]);
cb = colorbar(ax_m_4,'units','inches','Location', 'northoutside', 'position',fig_cfg.cb_size);
set(cb,'fontsize',15,'fontname','freeserif','TickDir','both');
cb_title=title(cb,'$$ \tau $$','fontsize', 22, 'Position', [1240, 0, 0]); % hor, ver, ?
set(cb_title, 'interpreter', 'latex');


%% plot using surf (Individual)
ax_m_5=subplot(7,7,5);
fig_cfg.ax_size = [loc_column_first + 12, loc_row_first, 5, 5];
period=1:5000;
fig_weig=1;
% ax_m_5=axes('fontsize',14, 'fontname','freeserif'); 
set(ax_m_5, 'Parent', fig_h)
set(ax_m_5,'Units','inches','Position',fig_cfg.ax_size);
set(ax_m_5,'fontsize',12);
grid on
hold on
% plots.X=HCST_traj.X.ens05.i0500(period);
plots.X=HCST_traj.X.i0500(8,period);
plots.Y=HCST_traj.Y.i0500(8,period);
plots.Z=HCST_traj.Z.i0500(8,period);
tmp.parula=parula(200+(max(period))); % from light blue, until darker yellow
plots=rmfield(plots,'C0');
for ti=1:length(plots.X)
    plots.C0(ti,1:2,1)=tmp.parula(ti+100,1);
    plots.C0(ti,1:2,2)=tmp.parula(ti+100,2);
    plots.C0(ti,1:2,3)=tmp.parula(ti+100,3);
end
 objsurf=surf(ax_m_5, [plots.X(:) plots.X(:)], [plots.Y(:) plots.Y(:)], [plots.Z(:) plots.Z(:)], ...
        plots.C0, ...  % Reshape and replicate data
     'FaceColor', 'none', ...    % Don't bother filling faces with color
     'EdgeColor', 'interp', ...  % Use interpolated color for edges
     'LineWidth', fig_weig, 'EdgeAlpha', 1);            % Make a thicker line
title1=title('(c) Individual', 'fontsize', 20);

view(210,23)
caxis([0 max(period)/100]);
xlabel('$$ x(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
xlim([-27 27])
ylabel('$$ y(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
ylim([-29 29]);
zlabel('$$ z(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
zlim([0 55]);
cb = colorbar(ax_m_5,'units','inches','Location', 'northoutside', 'position',fig_cfg.cb_size);
set(cb,'fontsize',15,'fontname','freeserif','TickDir','both');
cb_title=title(cb,'$$ \tau $$','fontsize', 22, 'Position', [1240, 0, 0]); % hor, ver, ?
set(cb_title, 'interpreter', 'latex');

%scatter
ax_m_6=subplot(7,7,6);
set(ax_m_6, 'Color','none');
set(ax_m_6, 'Parent', fig_h)
set(ax_m_6,'Units','inches','Position',fig_cfg.ax_size);
set(ax_m_6,'fontsize',12);
hold on

sc1=scatter3(ax_m_6,plots.X(1), plots.Y(1), plots.Z(1), 'MarkerFaceColor', plots.C0(1,1,:), ...
    'MarkerEdgeColor', 'none', 'LineWidth', 3);

view(210,23)
caxis([0 max(period)/100]);
xlabel('$$ x(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
xlim([-27 27])
ylabel('$$ y(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
ylim([-29 29]);
zlabel('$$ z(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
zlim([0 55]);
cb = colorbar(ax_m_6,'units','inches','Location', 'northoutside', 'position',fig_cfg.cb_size);
set(cb,'fontsize',15,'fontname','freeserif','TickDir','both');
cb_title=title(cb,'$$ \tau $$','fontsize', 22, 'Position', [1240, 0, 0]); % hor, ver, ?
set(cb_title, 'interpreter', 'latex');


% q3=norminv(0.75);
% q95=norminv(0.975);
% w95=(q95-q3)/(2*q3);

ax_m_7=subplot(7,7,7);
caxis([0 max(period)/100]);
fig_cfg.ax_size = [loc_column_first, loc_row_first-6, 17, 5];
set(ax_m_7, 'Parent', fig_h)
set(ax_m_7,'Units','inches','Position',fig_cfg.ax_size);

% boxplot(HCST.X.corr(:,1), 1, 'BoxStyle', 'outline', 'Colors', 'krrrbbbbbbbbbbbbbbbbbbbbb', 'ColorGroup', RCM_fig.tot_trend_rand_name, ...
%     'Whisker', w95, 'symbol','');

sc_1=scatter(ax_m_7, (1:cfg.pred_len), HCST.X.corr_em);
set(sc_1, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
set(sc_1, 'SizeData', 50.*fig_weig)

hold on
lineplot_1=plot(ax_m_7, (1:cfg.pred_len), HCST.X.corr_em, 'k-', 'linewidth',2);
set(get(get(lineplot_1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); %legend off

bp1=boxplot(ax_m_7, HCST.X.corr(:,1:cfg.pred_len), (1:cfg.pred_len).*0.1-0.1, 'BoxStyle', 'outline', 'symbol', '', 'Whisker', 2);
% set(bp1, 'Widths', 0.4)
set(bp1, {'linew'}, {1.5})

lineplot_2=plot(ax_m_7, (1:cfg.pred_len), HCST.X.corr_med, 'r-', 'linewidth',2);


%% brier score line
lineplot_3=plot(ax_m_7, (1:cfg.pred_len), bsm, 'magenta', 'linewidth',2);


% bp1=boxchart(double(HCST.X.corr(:,1:11)),  'symbol', '');

ylim([-0.2 1])
title1=title('(d) ACC(\tau)', 'fontsize', 20);
grid on

ylabel('ACC for $$ x(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
xlabel('lead $$ \tau $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)

set(ax_m_7,'fontsize',15);
set(ax_m_7,'YAxisLocation', 'right');

% box_vars = findall(fig_h,'Tag','Box');
% hLegend = legend(box_vars([3,2,4]), {'Group 1','Group 2','Group 3'});



%% 90% significance low & high range
conf = 0.90;
n=cfg.t_ini_max-cfg.t_ini_min+1;
alpha = 1 - conf;
pLo = alpha/2;
pUp = 1 - alpha/2;
crit = tinv([pLo pUp], n-1);
% xbar = mean(r); % = 0
xbar = 0; % = 0
r_crit=sqrt((crit.^2)./(n-2+(crit).^2));

line_sig=line(ax_m_7, 0:cfg.pred_len+1, repmat(r_crit(1), 1,cfg.pred_len+2), 'color', 'g', 'LineStyle', '--', 'linewidth',2);

%% 0.5, for Brier score
line_sig2=line(ax_m_7, 0:cfg.pred_len+1, repmat(0.5, 1,cfg.pred_len+2), 'color', [0.9290 0.6940 0.1250], 'LineStyle', '--', 'linewidth',2);



% [hLg]=legend(ax_m_7, [sc1, lineplot_2, line_sig], ...
%     {'ensmean', 'median(individual)', '90% sig.'}, ...
%     'Fontsize', 15.*fig_weig);  %% for median

[hLg]=legend(ax_m_7, [sc1, lineplot_2, lineplot_3, line_sig, line_sig2], ...
    {'ensmean', 'median(individual)', 'Brier score', '90% sig. (ACC)', '0.5'}, ...
    'Fontsize', 15.*fig_weig);  %% for median


% fig_cfg.ax_size = [loc_column_first, loc_row_first, 5, 5];
% set(ax_m_1,'Units','inches','Position',fig_cfg.ax_size);
% set(ax_m_2,'Units','inches','Position',fig_cfg.ax_size);
% fig_cfg.ax_size = [loc_column_first + 6, loc_row_first, 5, 5];
% set(ax_m_3,'Units','inches','Position',fig_cfg.ax_size);
% set(ax_m_4,'Units','inches','Position',fig_cfg.ax_size);
% fig_cfg.ax_size = [loc_column_first + 12, loc_row_first, 5, 5];
% set(ax_m_5,'Units','inches','Position',fig_cfg.ax_size);
% set(ax_m_6,'Units','inches','Position',fig_cfg.ax_size);
% fig_cfg.ax_size = [loc_column_first, loc_row_first-6, 17, 5];
% set(ax_m_7,'Units','inches','Position',fig_cfg.ax_size);

print(gcf, ['lorenz_trajectory_',num2str(ext_flag), '.png'], '-dpng');
















