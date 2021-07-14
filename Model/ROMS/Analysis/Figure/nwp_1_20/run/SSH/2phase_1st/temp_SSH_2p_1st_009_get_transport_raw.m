close all; clear all;  clc;
warning off;


% error('Error in summation of the KS, TS, SS');

% all_region2 ={'NWP'}
    close all;
    clearvars '*' -except regionind2 testnameind2 all_region2 all_testname2 all_var2
    % % % 
    % % % Read Model SST
    % % % interp
    % % % get RMS
    % % % get BIAS
    system_name=computer;
    if (strcmp(system_name,'PCWIN64'))
        % % for windows
        dropboxpath='C:\Users\user\Dropbox';
        addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
        addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
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

    % for snu_desktopd
%     testname=all_testname2{testnameind2}    % % need to change
    inputyear = [1985:2014];
%     inputyear1 = [2006:2015]; % % put year which you want to plot [year year ...]
%     inputyear2 = [2091:2100]; % % put year which you want to plot [year year ...]
    inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
    scenname ='historical';

    if (strcmp(system_name,'PCWIN64'))
        % % for windows
        cmip6dir = strcat('D:\Data\Model\CMIP6\NWP\'); % % where data files are
    elseif (strcmp(system_name,'GLNXA64'))
        cmip6dir ='/data1/kimyy/Model/CMIP6/';
    end

% translev = [0 4.5];
    
% start-------------------- get korea strait transport (CNRM-ESM2)
testname='CNRM-ESM2-1';
outputdir=[cmip6dir, 'transport',filesep, scenname, filesep, 'Omon', filesep, testname, filesep];
pngname = [outputdir, 'ES_transp_',num2str(inputyear(1)),'_',num2str(inputyear(end)),'.png'];
% if (exist(pngname , 'file') ~= 2)       
    if (exist(strcat(outputdir) , 'dir') ~= 7)
        mkdir(strcat(outputdir));
    end 

    for yearij=1:length(inputyear)
        tempyear=inputyear(yearij);
        yearstr=num2str(tempyear, '%04i');
	txtname=[outputdir, testname, '_monthly_', scenname, '_', yearstr, '.txt'];
	fid = fopen(txtname, 'w+');
	fprintf(fid, '%%Korea  	%%Tsugaru  %%soya \n');
	fclose(fid);
        for monthij=1:length(inputmonth)
            tempmonth=inputmonth(monthij);
            monthstr=num2str(tempmonth, '%02i');
            xData((12*(yearij-1))+monthij) = datenum([num2str(tempmonth,'%02i'),'-01-',num2str(tempyear,'%04i')]);

            vname='vo';
            uname='uo';
            vfiledir = strcat(cmip6dir, vname, filesep, scenname, filesep, 'Omon', filesep, testname, filesep); % % where data files are
            vfilename=[vfiledir, filesep, vname, '_Omon_', scenname, '_', testname, '_', yearstr, '.nc'];
            ufiledir = strcat(cmip6dir, uname, filesep, scenname, filesep, 'Omon', filesep, testname, filesep); % % where data files are
            ufilename=[ufiledir, filesep, uname, '_Omon_', scenname, '_', testname, '_', yearstr, '.nc'];
            lonfilename = [vfiledir, filesep, 'lon', '_Omon_', scenname, '_', testname, '.nc'];
            latfilename = [vfiledir, filesep, 'lat', '_Omon_', scenname, '_', testname, '.nc'];
            depthfilename = [vfiledir, filesep, 'depth', '_Omon_', scenname, '_', testname, '.nc'];
            tind=(tempmonth);
            
            ks_points = [20, 21, 28, 28]; % x1, x2, y1, y2
            ts_points = [31, 31, 37, 37]; % x1, x2, y1, y2
            ss_points = [33, 33, 43, 44]; % x1, x2, y1, y2

            v_ks = ncread(vfilename,vname,[ks_points(1) ks_points(3) 1 tind], ...
                [ks_points(2)-ks_points(1)+1 ks_points(4)-ks_points(3)+1 inf 1]);
            if (strcmp(testname, 'CNRM-ESM2-1')==1)
                v_ks = ncread(vfilename,vname,[ks_points(1) ks_points(3) 1 tind], ...
                [ks_points(2)-ks_points(1)+1 ks_points(4)-ks_points(3)+2 inf 1]);
                v_ks=mean(v_ks,2);
            else
                v_ks = ncread(vfilename,vname,[ks_points(1) ks_points(3) 1 tind], ...
                [ks_points(2)-ks_points(1)+1 ks_points(4)-ks_points(3)+1 inf 1]);
            end
            u_ts = ncread(ufilename,uname,[ts_points(1) ts_points(3) 1 tind], ...
                [ts_points(2)-ts_points(1)+1 ts_points(4)-ts_points(3)+1 inf 1]);  
            u_ss = ncread(ufilename,uname,[ss_points(1) ss_points(3) 1 tind], ...
                [ss_points(2)-ss_points(1)+1 ss_points(4)-ss_points(3)+1 inf 1]); 
            zname ='lev';
            zboundname ='lev_bounds';
            lev=ncread(depthfilename,zname);
            lev_bounds=ncread(depthfilename,zboundname);
%             diff_lev=diff(lev);
%             for thickij=1:length(lev)
%                 if thickij==1
%                     thick_lev(thickij)=lev(thickij)+diff_lev(thickij)/2.0;
%                 elseif thickij==length(lev)
%                     thick_lev(thickij)=diff_lev(thickij-1);
%                 else
%                     thick_lev(thickij)=(diff_lev(thickij-1)+diff_lev(thickij))/2.0;
%                 end
%             end

            thick_lev=lev_bounds(2,:)-lev_bounds(1,:);
            
            thick_lev_ks=repmat(thick_lev, [size(v_ks,1) 1]);
            thick_lev_ts=repmat(thick_lev, [size(u_ts,2) 1]);
            thick_lev_ss=repmat(thick_lev, [size(u_ss,2) 1]);
            res_thick_lev_ks = reshape(thick_lev_ks, [size(v_ks)]);
            res_thick_lev_ts = reshape(thick_lev_ts, [size(u_ts)]);
            res_thick_lev_ss = reshape(thick_lev_ss, [size(u_ss)]);
            
            vall_ks=res_thick_lev_ks.*v_ks;
            uall_ts=res_thick_lev_ts.*u_ts;
            uall_ss=res_thick_lev_ss.*u_ss;
            vall_ks(isnan(vall_ks))=0;
            vall_ks=mean(vall_ks,1);
            uall_ts(isnan(uall_ts))=0;
            uall_ts=mean(uall_ts,2);
            uall_ss(isnan(uall_ss))=0;
            uall_ss=mean(uall_ss,2);
            
            vsum_ks=sum(vall_ks(:),'omitnan');
            usum_ts=sum(uall_ts(:),'omitnan');
            usum_ss=sum(uall_ss(:),'omitnan');
            if (exist('lonw_ks')==0)
                lonboundname= 'bounds_lon';  % vertices : sw, se, ne, nw
                latboundname= 'bounds_lat';
%                 Korea Strait points
                lonw_ks=ncread(lonfilename,lonboundname,[1, ks_points(1),ks_points(3)], [1 1 1]);
                lone_ks=ncread(lonfilename,lonboundname,[2, ks_points(2),ks_points(4)], [1 1 1]);
                latw_ks=ncread(latfilename,latboundname,[1, ks_points(1),ks_points(3)], [1 1 1]);
                late_ks=ncread(latfilename,latboundname,[2, ks_points(2),ks_points(4)], [1 1 1]);
                xdist_ks=m_lldist([lonw_ks,lone_ks], [latw_ks,late_ks]) * 1000.0; %% (m)
%                 Tsugaru Strait points
                lons_ts=ncread(lonfilename,lonboundname,[1, ts_points(1), ts_points(3)], [1 1 1]);
                lonn_ts=ncread(lonfilename,lonboundname,[4, ts_points(2), ts_points(4)], [1 1 1]);
                lats_ts=ncread(latfilename,latboundname,[1, ts_points(1), ts_points(3)], [1 1 1]);
                latn_ts=ncread(latfilename,latboundname,[4, ts_points(2), ts_points(4)], [1 1 1]);
                ydist_ts=m_lldist([lons_ts,lonn_ts], [lats_ts,latn_ts]) * 1000.0; %% (m)
%                 Soya Strait points
                lons_ss=ncread(lonfilename,lonboundname,[1, ss_points(1), ss_points(3)], [1 1 1]);
                lonn_ss=ncread(lonfilename,lonboundname,[4, ss_points(2), ss_points(4)], [1 1 1]);
                lats_ss=ncread(latfilename,latboundname,[1, ss_points(1), ss_points(3)], [1 1 1]);
                latn_ss=ncread(latfilename,latboundname,[4, ss_points(2), ss_points(4)], [1 1 1]);
                ydist_ss=m_lldist([lons_ss,lonn_ss], [lats_ss,latn_ss]) * 1000.0; %% (m)
            end
            korea_tr((yearij-1)*12+monthij)=vsum_ks*xdist_ks/1e6;
            tsugaru_tr((yearij-1)*12+monthij)=usum_ts*ydist_ts/1e6;
            soya_tr((yearij-1)*12+monthij)=usum_ss*ydist_ss/1e6;
            
            disp([testname, ' ', num2str(tempyear),'Y ', num2str(tempmonth), 'M']);
            fid=fopen(txtname, 'a+');
	    fprintf(fid, '%8.3f %8.3f %8.3f \n', ...
            korea_tr((yearij-1)*12+monthij), tsugaru_tr((yearij-1)*12+monthij), soya_tr((yearij-1)*12+monthij));
	    fclose(fid);
        end
    end
%     clear lonw_ks
% end-------------------- get korea strait transport (CNRM-ESM2)


% clearvars '*' -except inputyear inputyear1 inputyear2 inputmonth scenname cmip6dir translev
