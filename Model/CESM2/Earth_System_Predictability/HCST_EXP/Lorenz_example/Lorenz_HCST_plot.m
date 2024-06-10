close all; clear all; clc;

dir.lroot = ['~/Desktop/backup/Research/Postdoc/03_IBS/', ...
    '2022_predictability_assimilation_run/Lorenz_predictability_example/Model_src-v1'];

dir.hcst = [dir.lroot, '/', 'ensPRED'];
dir.hcst_lqini = [dir.lroot, '/', 'ensPRED_lqini']; % low quality initialized
dir.ctl = [dir.lroot, '/', 'CTL'];

%% read control run
for t_ini=51:100
    str.t_ini=num2str(t_ini, '%04i');
    dir.ctl_tmp=[dir.ctl, '/', 'out'];
    fname_ctl=[dir.ctl_tmp, '/', 'out.', str.t_ini, '.nc'];
    CTL_X.(['p',str.t_ini])=ncread(fname_ctl, 'X');
    CTL_Y.(['p',str.t_ini])=ncread(fname_ctl, 'Y');
    CTL_Z.(['p',str.t_ini])=ncread(fname_ctl, 'Z');
end

%% read ensemble predictions
for t_ini=55:5:85
    str.t_ini=num2str(t_ini, '%04i');
    for ensnum=1:20
        str_ens=num2str(ensnum, '%02i');
        dir.hcst_tmp=[dir.hcst, '/', str.t_ini, '/', 'ens',str_ens, '/', 'out'];
        dir.hcst_lqini_tmp=[dir.hcst, '/', str.t_ini, '/', 'ens',str_ens, '/', 'out'];
        for t_pred=1:5
            str.t_hcst=num2str(t_ini+t_pred-1, '%04i');
            fname_hcst=[dir.hcst_tmp, '/', 'out.', str.t_ini, '.nc'];
            fname_hcst_lqini=[dir.hcst_lqini_tmp, '/', 'out.', str.t_ini, '.nc'];
            HCST_X.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=ncread(fname_hcst, 'X');
            HCST_Y.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=ncread(fname_hcst, 'Y');
            HCST_Z.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=ncread(fname_hcst, 'Z');
            HCST_lqini_X.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=ncread(fname_hcst_lqini, 'X');
            HCST_lqini_Y.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=ncread(fname_hcst_lqini, 'Y');
            HCST_lqini_Z.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=ncread(fname_hcst_lqini, 'Z');
        end
        disp([str.t_ini, ', ens : ', str_ens]);
    end
end

hold on
for t_ini=51:100
    str.t_ini=num2str(t_ini, '%04i');
    t_plot=t_ini:0.01:t_ini+0.99;
    plot(t_plot, CTL_X.(['p',str.t_ini]), 'b');
end
hold off


