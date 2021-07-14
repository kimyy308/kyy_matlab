%%== plot zosto yearly map
clearvars; close all;
set(groot,'DefaultAxesFontName','FreeSerif','DefaultAxesFontSize',14,'defaultColorbarFontsize',14);
set(groot,'DefaultTextFontName','FreeSerif','DefaultTextFontSize',16);
flag_savepng = 1;

VarName    = 'zosto';
VarUnit    = 'm';
Tables = ["Omon","interp"];
Experiments = ["historical","rcp85","rcp45"];
Models = ["IPSL-CM5A-LR","IPSL-CM5A-MR","MPI-ESM-LR","NorESM1-M"];

if( flag_savepng ), f1=figure('visible','off');
else,               f1=figure; end
a1=axes; hold on; grid on; box on;
for Table = Tables
for Experiment = Experiments
for Model = Models
    if( contains(Experiment,'rcp') )
        yr_range = (2006:2100);   % for year loop, historical era
    else
        yr_range = (1976:2005);   % for year loop, historical era
    end
    time = (yr_range(1)+1/24:1/12:yr_range(end)+1)';
    i_time = ((yr_range(1)-1976)*12+1:(yr_range(end)-1975)*12)';
    fprintf( '== %s %s %s %s %04i-%04i ==\n', VarName, Table, Experiment, Model, yr_range(1), yr_range(end) );

    path_fig = fullfile('figure',VarName,Experiment,Table,Model);
    if( flag_savepng && ~exist(path_fig,'dir') ), mkdir(path_fig); end
    path_in = fullfile(VarName,Experiment);
    fname_in = sprintf('%s_%s_%s_%s_%04i-%04i.nc',VarName,Table,Model,Experiment,yr_range(1),yr_range(end));
    file_in = fullfile(path_in,fname_in);

    lon  = ncread( file_in, 'lon' );
    lat  = ncread( file_in, 'lat' );
    if( length(lon) == numel(lon) ), [lon,lat] = ndgrid(lon,lat);   end
    % time = ncread( file_in, 'time' );
    var  = ncread( file_in, VarName );

    file_customZ=fullfile('zostoga',Experiment,sprintf('%s_%s_%s_%s_%04i-%04i.mat','zostoga','Omon',Model,Experiment,yr_range(1),yr_range(end)));
    customZ = load(file_customZ);
    
%     refpolygon = akp4polygon;
%     mask_model = inpolygon(lon,lat,refpolygon(:,1),refpolygon(:,2));
%     var = var.*mask_model;
%     dA  = dA.*mask_model;

    fprintf( '  Draw yearly map..' ); lap_time = tic; lap_time_l = tic;
    var_y = NaN( [size(lon),length(yr_range)] );
    for ii = 1:length(yr_range)
        nchar(1) = fprintf('  %d..',yr_range(ii)); elapsed = toc(lap_time_l);
        nchar(2) = fprintf(' %.0f sec. (%.0f sec. left)', elapsed, elapsed*(length(yr_range)-ii+1)/ii );
        fname_fig = sprintf('%s_%s_%s_%s_%04i.png',VarName,Table,Model,Experiment,yr_range(ii));
        i_time = (1:12) + (ii-1)*12;
        var_y(:,:,ii) = mean( var(:,:,i_time), 3 );
        c_lim = [-0.2 0.2]+mean(customZ.zostoga(i_time));
        
        lon_lim = [115,164]; %lon : 115 ~164
        lat_lim = [15,52];   %lat : 15 ~ 52
        tmp_var = var_y(:,:,ii);
        tmp_var( lon<lon_lim(1) | lon>lon_lim(2) | lat<lat_lim(1) | lat>lat_lim(2) ) = NaN;
        pcolor(lon,lat,tmp_var); shading flat;
        xlim(lon_lim); ylim(lat_lim);

%         [lon2,i_lon] = sort(lon,1);
%         lat2 = NaN(size(lat));
%         tmp_var2 = NaN(size(tmp_var));
%         for jj = 1:size(lat,2)
%             lat2(:,jj) = lat(i_lon(:,jj),jj);
%             tmp_var2(:,jj) = tmp_var(i_lon(:,jj),jj);
%         end
%         pcolor(lon2,lat2,tmp_var2); shading flat;

        cb=colorbar; caxis(c_lim); title(cb,'(m)');
        text(0.02,0.98,sprintf('%s\n%d',Model,yr_range(ii)),'verticalalignment','top','horizontalalignment','left','units','normalized');
        if( flag_savepng ), print(f1,fullfile(path_fig,fname_fig),'-dpng','-r200');
        else,               drawnow; end
        hold off;
        fprintf(repmat('\b',1,sum(nchar)));
    end     % for ii in yr_range
    fprintf('%7.1f sec\n', toc(lap_time) );
end     % for Models
end     % for Experiments
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

