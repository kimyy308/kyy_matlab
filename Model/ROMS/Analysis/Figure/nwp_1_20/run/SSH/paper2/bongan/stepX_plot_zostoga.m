%%== plot zostoga time series
clearvars; close all;
set(groot,'DefaultAxesFontName','FreeSerif','DefaultAxesFontSize',14,'defaultColorbarFontsize',14);
set(groot,'DefaultTextFontName','FreeSerif','DefaultTextFontSize',16);

VarName    = 'zosto';
VarUnit    = 'm';
Models = ["IPSL-CM5A-LR","IPSL-CM5A-MR","MPI-ESM-LR","NorESM1-M"];
path_fig = fullfile('figure','zostoga');

tot_var = NaN( length(Models), length(1976:2100)*12 );
tot_Z   = NaN( length(Models), length(1976:2100)*12 );
tot_Z_c = NaN( length(Models), length(1976:2100)*12 );
tot_customZ = NaN( length(Models), length(1976:2100)*12 );
tot_time = (1976+1/24:1/12:2101)';

for Table = ["Omon","interp"]
for Scenario = ["rcp85","rcp45"]
    fprintf( '%s %s %s..', 'zostoga', Table, Scenario ); lap_time = tic;
for Experiment = ["historical",Scenario]
    if( contains(Experiment,'rcp') )
        yr_range = (2006:2100);   % for year loop, historical era
    else
        yr_range = (1976:2005);   % for year loop, historical era
    end
    time = (yr_range(1)+1/24:1/12:yr_range(end)+1)';
    i_time = ((yr_range(1)-1976)*12+1:(yr_range(end)-1975)*12)';

    %--- read zosto and then average
    for Model = Models

        path_in = fullfile(VarName,Experiment);
        fname_in = sprintf('%s_%s_%s_%s_%04i-%04i.nc', ...
                            VarName,Table,Model,Experiment,yr_range(1),yr_range(end));
        file_in = fullfile(path_in,fname_in);

        lon  = ncread( file_in, 'lon' );
        lat  = ncread( file_in, 'lat' );
        if( length(lon) == numel(lon) ), [lon,lat] = ndgrid(lon,lat);   end
        % time = ncread( file_in, 'time' );
        var  = ncread( file_in, VarName );

        if( contains( Table, 'Omon' ) )
            path_A = fullfile('/data1/CMIP/cmip5/historical_extHDD/CMIP5/areacello/piControl/fx',Model);
            file_A = dir(fullfile(path_A,'areacello*'));
            dA  = ncread(fullfile(file_A.folder,file_A.name),'areacello');
        else
            dA  = cosd(lat);
        end
        oceanA = nansum( ~any(isnan(var),3).*dA, 'all');

        avg_var = squeeze(nansum( var.*dA, [1,2] )/oceanA);
        tot_var(ismember(Models,Model),i_time) = avg_var;

    end

    %--- read CMIP5 zostoga
    for Model = Models
        if( contains(Experiment,'rcp') )
            path_Z = fullfile('/data1/CMIP/cmip5/rcp_extHDD/CMIP5/zostoga',Experiment,'Omon',Model);
        else
            path_Z = fullfile('/data1/CMIP/cmip5/historical_extHDD/CMIP5/zostoga',Experiment,'Omon',Model);
        end
        file_Z = dir(fullfile(path_Z,'zostoga*'));
        Z = ncread(fullfile(file_Z.folder,file_Z.name),'zostoga',1,length(time));
        Z_c = Z - Z(1);
        if( contains(Experiment,'rcp') )
            Z_c = Z_c + 0.031218;
        end

        tot_Z(ismember(Models,Model),i_time) = Z;
        tot_Z_c(ismember(Models,Model),i_time) = Z_c;
    end

    %--- read calculated zostoga from CMIP5 thetao & so
    for Model = Models
        file_customZ=fullfile('zostoga',Experiment,sprintf('%s_%s_%s_%s_%04i-%04i.mat','zostoga','Omon',Model,Experiment,yr_range(1),yr_range(end)));
        customZ = load(file_customZ);
        tot_customZ(ismember(Models,Model),i_time) = customZ.zostoga;
    end

end     % for Experiments

title1 = sprintf('%s %s avg(zosto,global)',Scenario,Table);
title2 = sprintf('%s %s zostoga(provided)',Scenario,Table);
title3 = sprintf('%s %s zostoga(offset)',Scenario,Table);
title4 = sprintf('%s %s zostoga(calculated)',Scenario,Table);
f1=figure('visible','off'); a1=axes; hold on; grid on; box on; ylabel(['(',VarUnit,')']);
f2=figure('visible','off'); a2=axes; hold on; grid on; box on; ylabel(['(',VarUnit,')']);
f3=figure('visible','off'); a3=axes; hold on; grid on; box on; ylabel(['(',VarUnit,')']);
f4=figure('visible','off'); a4=axes; hold on; grid on; box on; ylabel(['(',VarUnit,')']);
for Model = Models
    plot(tot_time,tot_var(ismember(Models,Model),:),'linewidth',1,'parent',a1);
    plot(tot_time,tot_Z(ismember(Models,Model),:),'linewidth',1,'parent',a2);
    plot(tot_time,tot_Z_c(ismember(Models,Model),:),'linewidth',1,'parent',a3);
    plot(tot_time,tot_customZ(ismember(Models,Model),:),'linewidth',1,'parent',a4);
end
legend(a1,Models,'location','northwest'); title(a1,title1);
legend(a2,Models,'location','northwest'); title(a2,title2);
legend(a3,Models,'location','northwest'); title(a3,title3);
legend(a4,Models,'location','northwest'); title(a4,title4);
ylim(a1,[-0.1 0.6]); xlim(a1,[1975 2101]);
ylim(a2,[-0.1 0.6]); xlim(a2,[1975 2101]);
ylim(a3,[-0.1 0.6]); xlim(a3,[1975 2101]);
ylim(a4,[-0.1 0.6]); xlim(a4,[1975 2101]);
if( ~exist(path_fig,'dir') ), mkdir(path_fig); end
print(f1,fullfile(path_fig,replace(title1,' ','_')),'-dpng','-r200');
print(f2,fullfile(path_fig,replace(title2,' ','_')),'-dpng','-r200');
print(f3,fullfile(path_fig,replace(title3,' ','_')),'-dpng','-r200');
print(f4,fullfile(path_fig,replace(title4,' ','_')),'-dpng','-r200');

fprintf('%7.1f sec\n', toc(lap_time) );
end     % for Scenario
end     % for Tables

%% functions
function polygon = akp4polygon
% %  around Korean Peninsula polygon
    polygon = ...
           [117.0, 30.0;
            117.0, 52.0;
            142.5, 52;
            142.3, 47;
            142, 46.5;
            142, 45;
            142, 43;
            141, 43;
            141, 42.8;
            140.2, 42.6;
            140.2, 42.2;
            140.4, 41.8;
            140.5, 41;
            140.5, 38;
            137, 36;
            136, 35;
            133, 35;
            132, 34;
            131, 34;
            131, 33;
            131.0, 32.0;
            128.4, 30.0]; 
end

