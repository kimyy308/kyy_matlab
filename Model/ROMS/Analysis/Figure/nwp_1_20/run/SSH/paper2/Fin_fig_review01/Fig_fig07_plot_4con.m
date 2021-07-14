close all; clear all;  clc;
warning off;

all_testname1  = {'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'NorESM1-M', 'MPI-ESM-LR'};
all_testname2  = {'test53', 'test54', 'test55', 'test56'};
GCM_testnames = {'GCM-IPSL-L', 'GCM-IPSL-M', 'GCM-Nor', 'GCM-MPI'};
RCM_testnames = {'RCM-IPSL-L', 'RCM-IPSL-M', 'RCM-Nor', 'RCM-MPI'};
all_region2 ={'AKP4'};

variable = 'SSH';
% scenname='rcp45';
for regionind2=1:length(all_region2)
        close all;

        inputyear = [1993:2005]; % % put year which you want to plot [year year ...]
        yearstr_min=num2str(inputyear(1));
        yearstr_max=num2str(inputyear(end));
        inputmonth = [1:12]; % % put month which you want to plot [month month ...]
        system_name=computer;
        for folding=1:1
            if (strcmp(system_name,'PCWIN64'))
                % % for windows
                dropboxpath='C:\Users\User\Dropbox';
                addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
                addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
                addpath(genpath([dropboxpath '\source\matlab\Common\cptcmap']));
                addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
                addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
                addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Roms_tools\Preprocessing_tools']));
                addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
            elseif (strcmp(system_name,'GLNXA64'))
                dropboxpath='/home/kimyy/Dropbox';
                addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
                addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
                addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
                addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
            end
            
            run('nwp_polygon_point.m');
            regionname=all_region2{regionind2};
            switch(regionname)
            case('NWP')
                refpolygon=nwppolygon;
            case('AKP4') %% Around Korea Peninsula
                refpolygon=akp4polygon;
            end
            lonlat(1)=min(refpolygon(:,1));
            lonlat(2)=max(refpolygon(:,1));
            lonlat(3)=min(refpolygon(:,2));
            lonlat(4)=max(refpolygon(:,2));

            load('C:\Users\User\Dropbox\source\matlab\Common\Figure\gmt_ocean_mod2.mat')  % % set colormap (gmt_ocean, nonwhite)
        end
       
        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\paper2\fin_fig\'); % % where figure files will be saved
        elseif (strcmp(system_name,'GLNXA64'))
        end

        correction_right_fig=[-0.1600,0,0,0]; % right, up, width, height
        correction_upper_fig=[0,0,0,0];
        correction_large_fig=[0,0,0,0.020];
        hold on;
        tifname=strcat(figrawdir, 'fig07','.tif'); %% ~_year_month.jpg

%         f1=figure(1);
        byrmap = customcolormap_preset('red-yellow-blue');
        yrmap = byrmap(129:256,:);
        
% start-------------------- validiation point dot
        testname ='test53';
        filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
        tempyear=2005;
        param_script =['C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m'];
        run(param_script);
        opts = spreadsheetImportOptions("NumVariables", 12);
        opts.Sheet = "Sheet1";
        opts.DataRange = "B2:M20";
        opts.VariableNames = ["ofstation", "name", "lon", "lat", "M2amp", "S2amp", "K1amp", "O1amp", "M2phase", "S2phase", "K1phase", "O1phase"];
        opts.VariableTypes = ["double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
        opts = setvaropts(opts, "name", "WhitespaceRule", "preserve");
        opts = setvaropts(opts, "name", "EmptyFieldRule", "auto");
        tcon_obs = readtable("Z:\내 드라이브\Data\Observation\Tide_report_choi\4con_data4_obs_and_khoa.xlsx", opts, "UseExcel", false);
        clear opts
        
        Obs_data = tcon_obs;
        tconname=[filedir, num2str(tempyear,'%04i'), '\', testname, '_tcon_val_zeta_', num2str(tempyear, '%04i'), '.mat'];
        load(tconname)

        amplev=([0 420]);
        sb{1}=subplot(2,4,[1 2 5 6]);  % Auto-fitted to the figure.
        pos_sb{1}=get(sb{1}, 'pos'); % Get the position.delete(sb1); % Delete the subplot axes
        delete(sb{1}); % Delete the subplot axes
        ax{1}=axes;
        set(ax{1},'pos',pos_sb{1});
        m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
        m_gshhs_i('color',m_gshhs_line_color);
        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        hold on;
        amp_4con=sum(table2array(tcon_obs(:,5:8)),2);
        for linei=1:size(tcon_obs,1)
            if (amp_4con(linei)<=min(amplev))
                colind=1;
            elseif (amp_4con(linei)>=max(amplev))
                colind=128;
            else
                colind=round((amp_4con(linei)-min(amplev))/diff(amplev) *128.0);
            end
            m_line(table2array(tcon_obs(linei,4)),table2array(tcon_obs(linei,3)),'marker','o','color',yrmap(colind,:),'linewi',2,...
              'linest','none','markersize',8,'markerfacecolor',yrmap(colind,:), 'parent',ax{1});
        end
%                     m_plot(table2array(tcon_obs(:,4)),table2array(tcon_obs(:,3)),'marker','o','color','k', ...
%                     'markerfacecolor',bwrmap(colind,:)); 
        h{1} = colorbar;
        colormap(ax{1},yrmap);
        set(h{1},'fontsize',colorbar_fontsize);
        title(h{1},'(cm)','fontsize',colorbar_title_fontsize);
        caxis(ax{1},amplev);

        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  'box', m_grid_box_type, ...
                    'xticklabels', [120, 130, 140], 'xtick',[120, 130, 140], ...
                    'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'parent', ax{1});  
                
        txt_1_1=m_text(118, 51, '(A)', 'FontSize', m_grid_fontsize+5); 

        titlename = strcat('Obs station');  
        
        

                    
        opts = spreadsheetImportOptions("NumVariables", 12);
        opts.Sheet = "Sheet1";
        opts.DataRange = "B2:M20";
        opts.VariableNames = ["ofstation", "name", "lon", "lat", "M2amp", "S2amp", "K1amp", "O1amp", "M2phase", "S2phase", "K1phase", "O1phase"];
        opts.VariableTypes = ["double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
        opts = setvaropts(opts, "name", "WhitespaceRule", "preserve");
        opts = setvaropts(opts, "name", "EmptyFieldRule", "auto");
        tcon_obs = readtable("Z:\내 드라이브\Data\Observation\Tide_report_choi\4con_data4_obs_and_khoa.xlsx", opts, "UseExcel", false);
        clear opts

        for temp_testind=1:length(all_testname2)
            temp_testname=all_testname2{temp_testind};
            temp_filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\', temp_testname, '\run\'); % % where data files are
%                         temp_filename=[temp_filedir, num2str(tempyear,'%04i'), '\', temp_testname, '_harmonic_analysis_', varname, '_', regionname, '_', num2str(tempyear, '%04i'), '.nc'];
%                         temp_tcon=ncread(temp_filename, 'tcon');
%                         all_m2_amp(1:size(temp_tcon,1), 1:size(temp_tcon,2), temp_testind)=temp_tcon(:,:,tide_info.index(1),1).*100;
            temp_tconname=[temp_filedir, num2str(tempyear,'%04i'), '\', temp_testname, '_tcon_val_zeta_ref_', num2str(tempyear, '%04i'), '.mat'];
            load(temp_tconname)
            size1=size(tcon_model,1);
            size2=size(tcon_model,2);
            for cellx=1:size1
                for celly=1:4
                    all_tcon_model(size1*(temp_testind-1)+cellx,celly)=tcon_model{cellx,celly+4};
                    all_tcon_obs(size1*(temp_testind-1)+cellx,celly)=table2array(tcon_obs(cellx,celly+4));
                end
            end
        end

%                 if (exist('tcon_model') ~= 1)
        
        RCM_data = all_tcon_model(1:19,:);
        model_4con= sum(all_tcon_model,2);
        obs_4con =  sum(all_tcon_obs,2);
        
        
        sb{2}=subplot(2,5,[4 5 9 10]);  % Auto-fitted to the figure.
        pos_sb{2}=get(sb{2}, 'pos'); % Get the position.delete(sb1); 
        delete(sb{2}); % Delete the subplot axes
        ax{2}=axes;
        set(ax{2},'pos',pos_sb{2});
        
%                 num_sta=size(tcon_obs,1);
        mslplot{1}=scatter(obs_4con(1:cellx), model_4con(1:cellx), 'o', 'parent',ax{2});
        hold on
        mslplot{2}=scatter(obs_4con(cellx+1:cellx*2), model_4con(cellx+1:cellx*2), '+', 'parent',ax{2});
        mslplot{3}=scatter(obs_4con(cellx*2+1:cellx*3), model_4con(cellx*2+1:cellx*3), '*', 'parent',ax{2});
        mslplot{4}=scatter(obs_4con(cellx*3+1:cellx*4), model_4con(cellx*3+1:cellx*4), 'x', 'parent',ax{2});
        maxval=max([obs_4con,  model_4con]);
        mslplot{5}=plot(0:maxval, 0:maxval,'r', 'parent',ax{2});
        hold off
        
        for corri=1:4
            [corrcf, pval] = corrcoef(obs_4con(19*(corri-1)+1 : 19*(corri))', model_4con(19*(corri-1)+1 : 19*(corri))');
            tempcf(corri)=corrcf(1,2);
            temppval(corri)=pval(1,2);
        end
        
        xlabel(ax{2},'Observed Tidal amplitude  (cm)')
        ylabel(ax{2},'RCM Tidal amplitude (cm)')
%                     title(['all test', ', 4con amp, ', num2str(tempyear,'%04i')])

        % datetick('x','yyyy','keepticks')
        axis tight;
        % ylim(meanplotlev2)
        set(mslplot{1},'LineWidth',2);
        set(mslplot{2},'LineWidth',2);
        set(mslplot{3},'LineWidth',2);
        set(mslplot{4},'LineWidth',2);
        set(mslplot{5},'LineWidth',5);
        set(mslplot{1},'SizeData',100);
        set(mslplot{2},'SizeData',100);
        set(mslplot{3},'SizeData',100);
        set(mslplot{4},'SizeData',100);

        set(gca,'FontSize',m_grid_fontsize);
        grid on
        
        txt_1_2=text(20, 370, '(B)', 'FontSize', m_grid_fontsize+5); 

%                     lgd=legend(all_testname2);
        lgd=legend({'RCM-IPSL-L', 'RCM-IPSL-M', 'RCM-Nor', 'RCM-MPI'});

        set(lgd,'FontSize',m_grid_fontsize);
        set(lgd,'Position',[0.75 0.20, 0.20, 0.03]);
        set(lgd,'Orientation','vertical');
        
        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width*2 paper_position_height*1.2]) 

%         saveas(gcf,tifname,'tif');
        print('-dtiff','-r500',tifname)

        hold off;
        close all;
        save ('Z:\내 드라이브\research\Ph_D_course\2020_SSH_CMIP5_decadal_variation_around_korea_using_downscaling_ocean_model\data_dryad\Data07_tidal_constitute.mat', ...
            'Obs_data', 'RCM_data', ...
            'GCM_testnames', 'RCM_testnames')
end