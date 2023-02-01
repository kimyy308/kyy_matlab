% %  Created 29-Jan-2023 by Yong-Yub Kim  
clc; clear all; close all;

%% set path
tmp.dropboxpath = '/Volumes/kyy_raid/kimyy/Dropbox';
tmp.fs=filesep;
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);

dirs.root='/Volumes/kyy_raid/kimyy/Model/CMIP/CMIP6_analysis/tos/NA_reg/feb/runmean_9';
dirs.figroot='/Volumes/kyy_raid/kimyy/Figure/CMIP/CMIP6_analysis/tos/NA_reg/feb/runmean_9';
dir_res=dir(dirs.root);

tmp.varn='tos';

tmp.fs=filesep;
tmp.f_num=length(dir_res);

tmp.time = 1954:2096;

for ind_f=1:tmp.f_num
    tmp.f_name=[dirs.root, tmp.fs, dir_res(ind_f).name];
    if strcmp(tmp.f_name(end-2:end), '.nc')==1
        tmp.spn=split(dir_res(ind_f).name, '_');
        tmp.data= squeeze(ncread(tmp.f_name, tmp.varn));
        tmp.data_det = Func_0028_detrend_linear_1d(tmp.data)';
        CMIP6_hist{ind_f}.model=tmp.spn{end-1};
        CMIP6_hist{ind_f}.data=tmp.data(8:57);
        CMIP6_hist{ind_f}.data_det=Func_0028_detrend_linear_1d(CMIP6_hist{ind_f}.data)';
%         CMIP6_hist{ind_f}.data_det=tmp.data_det(8:57);
        
%         strcmp(tmp.f_name(end-2:end), '.nc')==1
%             disp(CMIP6_hist{ind_f}.model)
        [CMIP6_hist{ind_f}.psd, CMIP6_hist{ind_f}.psd_sig, freq, CMIP6_hist{ind_f}.var] = ...
            Func_0027_get_PSD_siglev(CMIP6_hist{ind_f}.data_det);

        CMIP6_fut{ind_f}.model=tmp.spn{end-1};
        CMIP6_fut{ind_f}.data=tmp.data(94:143);
        CMIP6_fut{ind_f}.data_det=Func_0028_detrend_linear_1d(CMIP6_fut{ind_f}.data)';
%         CMIP6_fut{ind_f}.data_det=tmp.data_det(94:143);
        [CMIP6_fut{ind_f}.psd, CMIP6_fut{ind_f}.psd_sig, freq, CMIP6_fut{ind_f}.var] = ...
            Func_0027_get_PSD_siglev(CMIP6_fut{ind_f}.data_det);

        CMIP6_var_diff(ind_f,:)=CMIP6_fut{ind_f}.var-CMIP6_hist{ind_f}.var;
        
% % %         %% psd plot (hist)
% % %         plot(1./freq,CMIP6_hist{ind_f}.psd, '-o', 'linewidth', 2);
% % %         hold on
% % %         plot(1./freq,CMIP6_hist{ind_f}.psd_sig, '--*', 'linewidth', 2);
% % %         xlabel('Period(Y)');
% % %         ylabel('PSD');
% % %         title(CMIP6_hist{ind_f}.model)
% % %         legend('hist', 'sig-95%')
% % %         
% % %         set(gca, 'fontsize', 20);
% % %         grid minor
% % % 
% % %         dirs.figdir=[dirs.figroot, tmp.fs, 'psd'];
% % %         system(['mkdir -p ', dirs.figdir]);
% % %         tmp.figname=[dirs.figdir, tmp.fs, ...
% % %                 'CMIP6_hist_', CMIP6_hist{ind_f}.model, '_', 'psd_', tmp.varn, '.tif'];
% % %         set(gcf, 'PaperPosition', [0, 0, 6, 4]);
% % %         saveas(gcf,tmp.figname,'tif');
% % %         RemoveWhiteSpace([], 'file', tmp.figname);
% % %         close all;
% % % 
% % %         %% psd plot (fut)
% % %         plot(1./freq,CMIP6_fut{ind_f}.psd, '-o', 'linewidth', 2);
% % %         hold on
% % %         plot(1./freq,CMIP6_fut{ind_f}.psd_sig, '--*', 'linewidth', 2);
% % %         xlabel('Period(Y)');
% % %         ylabel('PSD');
% % %         title(CMIP6_fut{ind_f}.model)
% % %         legend('fut', 'sig-95%')
% % %         
% % %         set(gca, 'fontsize', 20);
% % %         grid minor
% % % 
% % %         dirs.figdir=[dirs.figroot, tmp.fs, 'psd'];
% % %         system(['mkdir -p ', dirs.figdir]);
% % %         tmp.figname=[dirs.figdir, tmp.fs, ...
% % %                 'CMIP6_fut_', CMIP6_hist{ind_f}.model, '_', 'psd_', tmp.varn, '.tif'];
% % %         set(gcf, 'PaperPosition', [0, 0, 6, 4]);
% % %         saveas(gcf,tmp.figname,'tif');
% % %         RemoveWhiteSpace([], 'file', tmp.figname);
% % %         close all;
% % % 
% % % 
% % % % %         %% var plot (fut)
% % % % %         plot(1./freq,CMIP6_fut{ind_f}.var, '-o', 'linewidth', 2);
% % % % %         hold on
% % % % %         plot(1./freq,CMIP6_fut{ind_f}.psd_sig, '--*', 'linewidth', 2);
% % % % %         xlabel('Period(Y)');
% % % % %         ylabel('Variance');
% % % % %         title(CMIP6_fut{ind_f}.model)
% % % % %         legend('fut', 'sig-95%')
% % % % %         
% % % % %         set(gca, 'fontsize', 20);
% % % % %         grid minor
% % % % % 
% % % % %         dirs.figdir=[dirs.figroot, tmp.fs, 'var'];
% % % % %         system(['mkdir -p ', dirs.figdir]);
% % % % %         tmp.figname=[dirs.figdir, tmp.fs, ...
% % % % %                 'CMIP6_fut_', CMIP6_hist{ind_f}.model, '_', 'var_', tmp.varn, '.tif'];
% % % % %         set(gcf, 'PaperPosition', [0, 0, 6, 4]);
% % % % %         saveas(gcf,tmp.figname,'tif');
% % % % %         RemoveWhiteSpace([], 'file', tmp.figname);
% % % % %         close all;
% % % 
% % % 
% % %         %% raw-det plot
% % %         plot(1:50,CMIP6_hist{ind_f}.data_det, 'linewidth', 2);
% % %         hold on
% % %         plot(1:50,CMIP6_fut{ind_f}.data_det, 'linewidth', 2);
% % %         xlabel('Year');
% % %         ylabel(['Temp(^oC)']);
% % %         title(CMIP6_hist{ind_f}.model)
% % %         legend('hist', 'fut')
% % %        
% % %         set(gca, 'fontsize', 20);
% % %         grid minor
% % % 
% % %         dirs.figdir=[dirs.figroot, tmp.fs, 'raw_det'];
% % %         system(['mkdir -p ', dirs.figdir]);
% % %         tmp.figname=[dirs.figdir, tmp.fs, ...
% % %                 'CMIP6_', CMIP6_hist{ind_f}.model, '_', 'raw_det', tmp.varn, '.tif'];
% % %         set(gcf, 'PaperPosition', [0, 0, 6, 4]);
% % %         saveas(gcf,tmp.figname,'tif');
% % %         RemoveWhiteSpace([], 'file', tmp.figname);
% % %         close all;
% % % 
% % %         %% raw plot
% % %         plot(1:50,CMIP6_hist{ind_f}.data, 'linewidth', 2);
% % %         hold on
% % %         plot(1:50,CMIP6_fut{ind_f}.data, 'linewidth', 2);
% % %         xlabel('Year');
% % %         ylabel(['Temp(^oC)']);
% % %         title(CMIP6_hist{ind_f}.model)
% % %         legend('hist', 'fut')
% % %        
% % %         set(gca, 'fontsize', 20);
% % %         grid minor
% % % 
% % %         dirs.figdir=[dirs.figroot, tmp.fs, 'raw'];
% % %         system(['mkdir -p ', dirs.figdir]);
% % %         tmp.figname=[dirs.figdir, tmp.fs, ...
% % %                 'CMIP6_', CMIP6_hist{ind_f}.model, '_', 'raw_', tmp.varn, '.tif'];
% % %         set(gcf, 'PaperPosition', [0, 0, 6, 4]);
% % %         saveas(gcf,tmp.figname,'tif');
% % %         RemoveWhiteSpace([], 'file', tmp.figname);
% % %         close all;


    end
end


% p=polyfit(tmp.time,tmp.data,1);
% tr1=p(1);
% data_fit=tr1*tmp.time + p(2);
% plot(tmp.data)
% hold on
% plot(data_fit)
% hold off

%% raw data plot
for ind_f=1:tmp.f_num
    tmp.f_name=[dirs.root, tmp.fs, dir_res(ind_f).name];
    if strcmp(tmp.f_name(end-2:end), '.nc')==1
        tmp.data= squeeze(ncread(tmp.f_name, tmp.varn));
        plot(tmp.time,tmp.data, 'linewidth', 2);
        hold on
        xlabel('Year');
        ylabel(['Temp(^oC)']);
        
        set(gca, 'fontsize', 20);
        grid minor
    end
end
%% specify figdir, figname, save
dirs.figdir=[dirs.figroot, tmp.fs, 'raw'];
system(['mkdir -p ', dirs.figdir]);
tmp.figname=[dirs.figdir, tmp.fs, ...
        'CMIP6_all', '_', 'raw_', tmp.varn, '.tif'];
set(gcf, 'PaperPosition', [0, 0, 8, 4]);
saveas(gcf,tmp.figname,'tif');
RemoveWhiteSpace([], 'file', tmp.figname);
close all;


%% detrended data plot
for ind_f=1:tmp.f_num
    tmp.f_name=[dirs.root, tmp.fs, dir_res(ind_f).name];
    if strcmp(tmp.f_name(end-2:end), '.nc')==1
        tmp.data= squeeze(ncread(tmp.f_name, tmp.varn));
        tmp.data_det = Func_0028_detrend_linear_1d(tmp.data)';
        plot(tmp.time,tmp.data_det, 'linewidth', 2);
        hold on
        xlabel('Year');
        ylabel(['Temp(^oC)']);
        
        set(gca, 'fontsize', 20);
        grid minor
    end
end
%% specify figdir, figname, save
dirs.figdir=[dirs.figroot, tmp.fs, 'raw'];
system(['mkdir -p ', dirs.figdir]);
tmp.figname=[dirs.figdir, tmp.fs, ...
        'CMIP6_all', '_', 'raw_det_', tmp.varn, '.tif'];
set(gcf, 'PaperPosition', [0, 0, 8, 4]);
saveas(gcf,tmp.figname,'tif');
RemoveWhiteSpace([], 'file', tmp.figname);
close all;


%% spectrum plot (all_hist)
for ind_f=1:tmp.f_num
    tmp.f_name=[dirs.root, tmp.fs, dir_res(ind_f).name];
    if strcmp(tmp.f_name(end-2:end), '.nc')==1
       
        plot(1./freq,CMIP6_hist{ind_f}.psd, '-o', 'linewidth', 2);
        hold on
        xlabel('Period(Y)');
        ylabel('PSD');
        
        set(gca, 'fontsize', 20);
        grid minor
    end
end
axis tight
ylim([0 30])

%% specify figdir, figname, save
dirs.figdir=[dirs.figroot, tmp.fs, 'psd'];
system(['mkdir -p ', dirs.figdir]);
tmp.figname=[dirs.figdir, tmp.fs, ...
        'CMIP6_all', '_', 'hist_', tmp.varn, '.tif'];
set(gcf, 'PaperPosition', [0, 0, 6, 4]);
saveas(gcf,tmp.figname,'tif');
RemoveWhiteSpace([], 'file', tmp.figname);
close all;

%% spectrum plot (all_future)
for ind_f=1:tmp.f_num
    tmp.f_name=[dirs.root, tmp.fs, dir_res(ind_f).name];
    if strcmp(tmp.f_name(end-2:end), '.nc')==1
       
        plot(1./freq,CMIP6_fut{ind_f}.psd, '-o', 'linewidth', 2);
        hold on
        xlabel('Period(Y)');
        ylabel('PSD');
        
        set(gca, 'fontsize', 20);
        grid minor
    end
end
axis tight
ylim([0 30])

%% specify figdir, figname, save
dirs.figdir=[dirs.figroot, tmp.fs, 'psd'];
system(['mkdir -p ', dirs.figdir]);
tmp.figname=[dirs.figdir, tmp.fs, ...
        'CMIP6_all', '_', 'fut_', tmp.varn, '.tif'];
set(gcf, 'PaperPosition', [0, 0, 6, 4]);
saveas(gcf,tmp.figname,'tif');
RemoveWhiteSpace([], 'file', tmp.figname);
close all;




tmp.signum_hist(1:25)=0;
tmp.sigvar_hist(1:25)=0;
%% spectrum plot (all_hist), statistically significance
for ind_f=1:tmp.f_num
    tmp.f_name=[dirs.root, tmp.fs, dir_res(ind_f).name];
    if strcmp(tmp.f_name(end-2:end), '.nc')==1
        tmp.n=length(CMIP6_hist{ind_f}.psd);
        for i=1:tmp.n
            if (CMIP6_hist{ind_f}.psd(i) < CMIP6_hist{ind_f}.psd_sig(i))
                tmp.psd(i)=NaN;
            else
                tmp.psd(i)=CMIP6_hist{ind_f}.psd(i);
            end
            if (isfinite(tmp.psd(i)))
                tmp.signum_hist(i)=tmp.signum_hist(i)+1;
                tmp.sigvar_hist(i)=tmp.sigvar_hist(i)+tmp.psd(i).*freq(i);
            end
        end
        
        

        plot(1./freq,tmp.psd, '-o', 'linewidth', 2);
        hold on
        xlabel('Period(Y)');
        ylabel('PSD');
        
        set(gca, 'fontsize', 20);
        grid minor
    end
end
tmp.sigvar_hist=tmp.sigvar_hist./tmp.signum_hist

axis tight
ylim([0 30])

%% specify figdir, figname, save
dirs.figdir=[dirs.figroot, tmp.fs, 'psd'];
system(['mkdir -p ', dirs.figdir]);
tmp.figname=[dirs.figdir, tmp.fs, ...
        'CMIP6_all', '_sig_', 'hist_', tmp.varn, '.tif'];
set(gcf, 'PaperPosition', [0, 0, 6, 4]);
saveas(gcf,tmp.figname,'tif');
RemoveWhiteSpace([], 'file', tmp.figname);
close all;


tmp.signum_fut(1:25)=0;
tmp.sigvar_fut(1:25)=0;

%% spectrum plot (all_fut), statistically significance
for ind_f=1:tmp.f_num
    tmp.f_name=[dirs.root, tmp.fs, dir_res(ind_f).name];
    if strcmp(tmp.f_name(end-2:end), '.nc')==1
        tmp.n=length(CMIP6_fut{ind_f}.psd);
        for i=1:tmp.n
            if (CMIP6_fut{ind_f}.psd(i) < CMIP6_fut{ind_f}.psd_sig(i))
                tmp.psd(i)=NaN;
            else
                tmp.psd(i)=CMIP6_fut{ind_f}.psd(i);
            end
            if (isfinite(tmp.psd(i)))
                tmp.signum_fut(i)=tmp.signum_fut(i)+1;
                tmp.sigvar_fut(i)=tmp.sigvar_fut(i)+tmp.psd(i).*freq(i);                
            end
        end
        
        plot(1./freq,tmp.psd, '-o', 'linewidth', 2);
        hold on
        xlabel('Period(Y)');
        ylabel('PSD');
        
        set(gca, 'fontsize', 20);
        grid minor
    end
end
tmp.sigvar_fut=tmp.sigvar_fut./tmp.signum_fut

axis tight
ylim([0 30])

%% specify figdir, figname, save
dirs.figdir=[dirs.figroot, tmp.fs, 'psd'];
system(['mkdir -p ', dirs.figdir]);
tmp.figname=[dirs.figdir, tmp.fs, ...
        'CMIP6_all', '_sig_', 'fut_', tmp.varn, '.tif'];
set(gcf, 'PaperPosition', [0, 0, 6, 4]);
saveas(gcf,tmp.figname,'tif');
RemoveWhiteSpace([], 'file', tmp.figname);
close all;





% plot(CMIP6_var_diff(3:36,2))
% hold on
% line(1:35,repmat(0, [35,1]))
% hold off
% ylim([-0.5 0.5])