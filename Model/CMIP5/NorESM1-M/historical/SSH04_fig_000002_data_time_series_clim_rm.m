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
% var_names2 = {'zos'}; 
var_names2 = {'thetao'}; 

model_name = 'NorESM1-M';
scen_name = 'historical';

for nvar2 = 1:length(var_names2)
    clearvars '*' -except workdir nvar2 var_names2 model_name scen_name
    var_name = var_names2{nvar2};
    
    %path in the damo server
    cseof_output_dir = [workdir, '/cseofs/', model_name, '/', scen_name, '/', var_name]
    load([cseof_output_dir, '/lv_layer_merge.mat']);
    var_name = var_names2{nvar2};
    
    for nyear = 1:length(inputyear)
        tempyear = inputyear(nyear);
        for month=1:12
            xData((12*(nyear-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
        end
    end
    
    data3=reshape(data,size(data,1),12,size(data,2)/12);
    clim_data=squeeze(mean(data3,3));
    % plot(mean(data,1)
    clim_data_mean=mean(clim_data,1);
    clear clim_data2 clim_data3
    for y=1:30
        for i=1:12
            clim_data2(:,i,y)=data3(:,i,y)-clim_data(:,i);
        end
    end
    clim_data3=reshape(clim_data2,size(data,1),360);
    
    switch(var_name)
        case('zos')
            mvar=mean(clim_data3,1)*1000.0;
        otherwise
            mvar=mean(clim_data3,1);
    end
    p=polyfit(xData,mvar,1);
    p2=polyfit(1:length(mvar),mvar,1);
    mvar2=xData*p(1)+p(2);
    
    p=polyfit(xData,mvar,1);
    p2=polyfit(1:length(mvar),mvar,1);
    mvar2=xData*p(1)+p(2);

    mvarplot=plot(xData,mvar,'k')
    hold on
    mvarplot2=plot(xData,mvar2,'Color','r')
    
    figure_output_dir = [workdir, '/figure/', model_name, '/', scen_name, '/', var_name]
    tifname=[figure_output_dir, ...
        '/SSH_04_000002_NWP_tas_clim_removed_time_series_', ...
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
    ylabel(['Mean ', var_name, '(atm) (', units,')'])
    title(['NWP Mean ', var_name, '(', num2str(min(inputyear),'%04i'), ...
        '-',num2str(max(inputyear),'%04i'),'), ', ...
        num2str(round(p2(1)*12.0,3)), ' ', units,'/y'])
    datetick('x', 'yy', 'keepticks')
    axis tight;
    % ylim(meanplotlev)
    set(mvarplot, 'LineWidth', 2);
    set(gca, 'FontSize', 15);
    grid on
    grid minor
    hold off
    saveas(gcf, tifname, 'tif');
    close all;
end
