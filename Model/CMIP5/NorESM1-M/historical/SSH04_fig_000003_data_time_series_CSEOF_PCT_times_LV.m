clc;clear all;close all;

% This is the code to convert acsii file from binary file, 
%fit to cseof analysis
system_name=computer;
if (strcmp(system_name,'PCWIN64'))
    % % for windows
    dropboxpath='C:\Users\KYY\Dropbox';
    addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
    addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
elseif (strcmp(system_name,'GLNXA64'))
    dropboxpath='/home/kimyy/Dropbox';
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
end

presentdir = pwd;

workdir='/data1/kimyy/etc/CMIP5_CSEOF';
% var_names2 = {'tas', 'psl', 'hur', 'rsds', 'zos'};  % var, psl, hur, rsds, ua, va
% var_names2 = {'ua', 'va'};  
var_names2 = {'zos'};  % var, psl, hur, rsds, ua, va
% var_names2 = {'thetao'};

model_name = 'NorESM1-M';
scen_name = 'historical';

for nvar2 = 1:length(var_names2)
    clearvars '*' -except workdir nvar2 var_names2 model_name scen_name
    var_name = var_names2{nvar2};
    
    %path in the damo server
    cseof_output_dir = [workdir, '/cseofs/', model_name, '/', scen_name, '/', var_name]
    load([cseof_output_dir, '/lv_layer_merge.mat']);
    var_name = var_names2{nvar2};
    
    mode = 3;
    %load cseof pc time series
    pct_data = importdata(cpct_tt)';
    pct = pct_data(:)';
    pct = reshape(pct,tlen,length(pct)/tlen); %[time mode]

    LV_increment1(LV_increment1==0) = NaN;

    for nyear = 1:length(inputyear)
        tempyear = inputyear(nyear);
        for month=1:12
            xData((12*(nyear-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
        end
    end

    % % % % % %  LV_increment --> %% [layer x y period mode]
    mean_LV_increment1=squeeze(mean(mean(LV_increment1,2,'omitnan'),3,'omitnan'));

    for i = 1 : mode
        for cycle = 1:tlen/T  %% tlen/period
            for period = 1:T  %% nested period
                combined_pct_LV(T*(cycle-1)+period,i)=mean_LV_increment1(period,i)*pct(T*(cycle-1)+period,i);
            end
        end
        switch(var_name)
            case('zos')
                mtas=squeeze(combined_pct_LV(:,i))'.*1000.0;
            otherwise
                mtas=squeeze(combined_pct_LV(:,i))';
        end
        p=polyfit(xData,mtas,1);
        p2=polyfit(1:length(mtas),mtas,1);
        mtas2=xData*p(1)+p(2);

        mtasplot=plot(xData,mtas,'k')
        hold on
        mtasplot2=plot(xData,mtas2,'Color','r')
        figure_output_dir = [workdir, '/figure/', model_name, '/', scen_name, '/', var_name]
        tifname=[figure_output_dir, ...
            '/SSH_04_000003_NWP_PCT_times_LV_time_series_mode_',num2str(i,'%02i'),'_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.tif'];
        xlabel('year')
        switch(var_name)
            case('tas')
                units='^oC';
            case('psl')
                units ='pa';
            case('hur')
                units ='%';
            case('rsds')
                units ='w/m^2';
            case('zos')
                units ='mm';
            case('thetao')
                units='^oC';
            case('ua')
                units ='m/s';
            case('va')
                units ='m/s';
            otherwise
                ('?')
        end
        ylabel(['Mean ', var_name, ' (', units,')'])
%         title(['NWP Mean ', var_name, ' PCT*LV mode ', num2str(i,'%02i'), '(', num2str(min(inputyear),'%04i'), ...
%             '-',num2str(max(inputyear),'%04i'),'), ', ...
%             num2str(round(p2(1)*12.0,3)), ' ', units, '/y'])
         title(['', var_name, ' PCT*LV mode ', num2str(i,'%02i'), ', ', ...
            num2str(round(p2(1)*12.0, 3)), ' ', units, '/y'])
        datetick('x', 'yy', 'keepticks')
        axis tight;
        % ylim(meanplotlev)
        set(mtasplot, 'LineWidth', 2);
        set(gca, 'FontSize', 15);
        grid on
        grid minor
        hold off
        saveas(gcf, tifname, 'tif'); 
    end
end
