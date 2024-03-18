clearvars; close all;
S = shaperead('landareas.shp');
% name_rgn = 'Glob';  map_proj = 'robinson';  x_lim = [0 360];    y_lim = [-80 89];   fig_size = [0,0,6,2.8]; ax_size = [-0.15,-0.2,5.4,3];   cb_size = [5.15,0.20,0.15,2.3]; title_pos = [0.5,1.0];
% name_rgn = 'NP';    map_proj = 'robinson';  x_lim = [98 284];   y_lim = [-20 65];   fig_size = [0,0,6,2.8]; ax_size = [ 0.15,-0.1,5.1,2.6]; cb_size = [5.15,0.20,0.15,2.3]; title_pos = [0.5,1.05];
% name_rgn = 'NWP';   map_proj = 'robinson';  x_lim = [113 165];  y_lim = [15 58];    fig_size = [0,0,6,2.8]; ax_size = [-0.15,-0.2,5.4,3];   cb_size = [5.15,0.20,0.15,2.3]; title_pos = [0.5,1.05];
name_rgn = 'Glob';  map_proj = 'eqdcylin';  x_lim = [0 360];    y_lim = [-80 89];   fig_size = [0,0,6,2.8]; ax_size = [-0.15,0.15,5.4,2.7];	cb_size = [5.15,0.20,0.15,2.3]; title_pos = [0.5,0.93];
% name_rgn = 'NP';    map_proj = 'eqdcylin';  x_lim = [98 284];   y_lim = [-20 65];   fig_size = [0,0,6,2.8]; ax_size = [ 0.15,0.15,5.0,2.6]; cb_size = [5.15,0.20,0.15,2.3]; title_pos = [0.5,0.90];
% name_rgn = 'NWP';   map_proj = 'eqdcylin';  x_lim = [113 165];  y_lim = [15 58];    fig_size = [0,0,6,2.8]; ax_size = [ 0.05,0.05,5.4,2.4];	cb_size = [5.15,0.20,0.15,2.3]; title_pos = [0.5,1.06];

Varname = 'tos';    VarScale = 1;	Varunit = '\circC';     VarNC = Varname;                c_lim = [-2,33];        c_map = parula(255);
models_rewrap = readcell('rewrap_models.txt');
models_unwrap = readcell('unwrap_models.txt');
models_wrap   = readcell('wrap_models.txt');
files = dir('/Volumes/kyy_raid/kimyy/bongan/example/data/*.nc');

for i_model = 1:length(files)
% for i_model = [1:2]
    file_parts = split(files(i_model).name,{'_','.'});
    Table = file_parts{2};
    Model = file_parts{3};
    Grd   = file_parts{6};
    Time  = split(file_parts{end-1},'-');
    Time  = Time{1};
    
    fprintf('%02d_%s_%s  ',i_model,Model,Grd); lap_time = tic;

    file_in = fullfile(files(i_model).folder,files(i_model).name);
    finfo = ncinfo( file_in );
    lon_src = double(ncread_var( file_in, {'longitude','lon','nav_lon'} ));
    lat_src = double(ncread_var( file_in, {'latitude','lat','nav_lat'} ));
    VarDim = length(finfo.Variables(ismember({finfo.Variables.Name},Varname)).Dimensions);
    if VarDim < 3, continue;    end
    var_src = ncread(file_in,VarNC,[1,1,1],[inf,inf,1]);

    if( length( lon_src ) == numel( lon_src ) )     % if lon_src is a 1-D array
        [lon_src, lat_src] = ndgrid( lon_src, lat_src );
    end
    [n_lon_src, n_lat_src] = size( lon_src );
    
    X = lon_src; Y = lat_src; C = var_src;
    % need to do something
    if ~isempty(ismember(Model,models_rewrap))
%         ncks -O --msa -d $dmn,-2 -d $dmn,1,-2 $FILE_TMP  $FILE_TMP
        fprintf('rewrap  ');
        X = X([end-1,2:end-1],:); Y = Y([end-1,2:end-1],:); C = C([end-1,2:end-1],:);
    elseif ~isempty(ismember(Model,models_unwrap))
%         ncks -O -d $dmn,0,-2  $FILE_TMP  $FILE_TMP
        fprintf('unwrap  ');
        X = X(1:end-1,:); Y = Y(1:end-1,:); C = C(1:end-1,:);
    elseif ~isempty(ismember(Model,models_wrap))
%         ncks -O --msa -d $dmn,-1 -d $dmn,0,-1  $FILE_TMP  $FILE_TMP
        fprintf('wrap    ');
        X = X([end,1:end],:); Y = Y([end,1:end],:); C = C([end,1:end],:);
    else
        fprintf('org     ');
    end
    
    fig_name = sprintf('%02d_%s_%s_%s_%s',i_model,Model,Grd,Table,Time);
        fig_h = figure('name',fig_name,'PaperUnits','inches','PaperPosition',fig_size,'position',fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
        %-- map setting
        ax_m = axesm('MapProjection',map_proj,'grid','off','fontsize',12,'fontname','freeserif'); axis off; hold on;
        setm(ax_m,'origin',[0,mean(x_lim)],'MapLatLimit',y_lim,'MapLonLimit',x_lim);
        setm(ax_m,'ParallelLabel','off','MeridianLabel','off','PLabelMeridian','west','MLabelParallel','south');
        set(ax_m,'Units','inches','Position',ax_size);
        text(ax_m,title_pos(1),title_pos(2),fig_name,'units','normalized','horizontalalignment','center','verticalalignment','middle','fontsize',14,'fontname','freeserif','interpreter','none')
        %-- caxis & colorbar
        caxis(ax_m,c_lim/VarScale); colormap(c_map);
        cb = colorbar(ax_m,'units','inches','position',cb_size);
        set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
        ylabel(cb,[replace(VarNC,'_',' '), '(',Varunit,')']);
        %-- draw on ax_m
        h_pc = pcolorm(Y,X,C/VarScale,'parent',ax_m); shading flat;
        geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
        setm(ax_m,'frame','on','FLineWidth',1);

        fig_dir = [mfilename,'/',Varname,'_',name_rgn];
        if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
        print(fig_h, [fig_dir,'/',fig_name], '-dpng');

    fprintf('%7.1f sec\n', toc(lap_time) );

end     % for i_model

%% functions

function finfo = select_file( list, yr )

    finfo = struct([]);
    for ll = 1 : length( list )
        fname = list(ll).name;
        fname_split = strsplit( fname, {'_','.'} );
        fyear_str   = strsplit( fname_split{end-1}, '-' );
        fyear_start = str2num( fyear_str{1}(1:4) );
        if( length( fyear_str ) > 1 )
            fyear_end   = str2num( fyear_str{2}(1:4) );
        else
            fyear_end   = fyear_start;
        end
        if( yr >= fyear_start && yr <= fyear_end )
            finfo = struct;
            finfo.folder     = list(ll).folder;
            finfo.fname      = list(ll).name;
            finfo.Table      = fname_split{2};
            finfo.Model      = fname_split{3};
            finfo.Experiment = fname_split{4};
            finfo.Variant    = fname_split{5};
            finfo.Grid       = fname_split{6};
            finfo.fyear_start = fyear_start;
            finfo.fyear_end  = fyear_end;
            return;
        end
    end

end

function [var,varname] = ncread_var( fname, varnames, varargin)
    finfo   = ncinfo( fname );
    kk = find(ismember({finfo.Variables.Name},varnames),1,'last');
    if( ~isempty(kk) ),	varname = finfo.Variables(kk).Name;
    else,               error('Can''t find variable name');
    end
    if( nargin > 2 ),	var = ncread( fname, varname, varargin{:} );
    else,               var = ncread( fname, varname );
    end
end

function dir_var_labs = sort_variant_labels( dir_var_labs )
    tmp_str = struct([]);
    for i_var_labs = 1:length(dir_var_labs)
        tmp = squeeze(split(dir_var_labs(i_var_labs).name,{'r','i','p','f'}));
        tmp_str(i_var_labs).r = str2double(tmp{2});
        tmp_str(i_var_labs).i = str2double(tmp{3});
        tmp_str(i_var_labs).p = str2double(tmp{4});
        tmp_str(i_var_labs).f = str2double(tmp{5});
    end
    [~,i_sort] = sortrows(struct2table(tmp_str),{'r','i','p','f'});
    dir_var_labs = dir_var_labs(i_sort);
    return
end

function dir_grds = sort_grid_labels( dir_grds )
    dir_grds(contains({dir_grds.name},'_')) = [];
    return
end

