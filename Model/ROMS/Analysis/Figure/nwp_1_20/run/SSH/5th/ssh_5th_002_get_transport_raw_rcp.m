close all; clear all;  clc;
warning off;

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
        dropboxpath='C:\Users\KYY\Dropbox';
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
    inputyear = [2006:2100];
    inputyear1 = [2006:2015]; % % put year which you want to plot [year year ...]
    inputyear2 = [2091:2100]; % % put year which you want to plot [year year ...]
    inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
    scenname ='rcp45';

    if (strcmp(system_name,'PCWIN64'))
        % % for windows
        cmip5dir = strcat('G:\Data\Model\CMIP5\'); % % where data files are
    elseif (strcmp(system_name,'GLNXA64'))
        cmip5dir ='/data1/CMIP/cmip5/rcp_ocn_extHDD/CMIP5/'
    end

translev = [0 4.5];
    
% start-------------------- get korea strait transport (IPSL-CM5A-LR)
testname='IPSL-CM5A-LR';
outputdir=[cmip5dir, 'transport/', scenname, '/Omon/', testname, '/'];
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
	fprintf(fid, '%%Korea  	 \n');
	fclose(fid);
        for monthij=1:length(inputmonth)
            tempmonth=inputmonth(monthij);
            monthstr=num2str(tempmonth, '%02i');
            xData((12*(yearij-1))+monthij) = datenum([num2str(tempmonth,'%02i'),'-01-',num2str(tempyear,'%04i')]);

            varname='vo';
            filedir = strcat(cmip5dir, varname, '/', scenname, '/Omon/', testname, '/'); % % where data files are
            flag_file_in = false;
            list = dir( [ filedir, '/', varname, '*' ]); 
            for kk = 1 : length( list )
                fname_in    = list(kk).name;
                fname_split = strsplit( fname_in, {'_','.'} );
                fyear_str   = strsplit( fname_split{end-1}, '-' );
                fyear_start = str2num( fyear_str{1}(1:4) );
                fyear_end   = str2num( fyear_str{2}(1:4) );
                if( tempyear >= fyear_start && tempyear <= fyear_end &&     ...                 
                        strcmp( fname_split{2}, 'Omon' ) &&         ...
                        strcmp( fname_split{3}, testname ) &&      ...                 
                        strcmp( fname_split{4}, scenname ) )
                    flag_file_in = true;            break;
                end         
            end         
            if( ~flag_file_in )
                fprintf('Source File for %04i does not Exist. Continue...\n',year);   continue;
            end
            vfilename=[filedir, '/', fname_in];
            tind_v=(tempyear-fyear_start)*12+tempmonth;

            v = ncread(vfilename,varname,[27 100 1 tind_v], [1 1 inf 1]);  %% cut horizontal area [x,y,z] (wider than target area)
            lev=ncread(vfilename,'lev');
            diff_lev=diff(lev);
            for thickij=1:length(lev)
                if thickij==1
                    thick_lev(thickij)=lev(thickij)+diff_lev(thickij)/2.0;
                elseif thickij==length(lev)
                    thick_lev(thickij)=diff_lev(thickij-1);
                else
                    thick_lev(thickij)=(diff_lev(thickij-1)+diff_lev(thickij))/2.0;
                end
            end
            vsum=sum(squeeze(v).*thick_lev','omitnan');
            if (exist('lons')==0)
                lons=ncread(vfilename,'lon',[26,100], [1 1]);
                lone=ncread(vfilename,'lon',[28,100], [1 1]);
                lats=ncread(vfilename,'lat',[26,100], [1 1]);
                late=ncread(vfilename,'lat',[28,100], [1 1]);
                xdist=m_lldist([lons,lone], [lats,late]) /2.0 * 1000.0; %% (m)
            end
            korea_tr((yearij-1)*12+monthij)=vsum*xdist/1000000.0;
            disp([testname, num2str(tempyear),'Y', num2str(tempmonth), 'M']);
            fid=fopen(txtname, 'a+');
	    fprintf(fid, '%8.3f \n', korea_tr((yearij-1)*12+monthij));
	    fclose(fid);
        end
    end
    % % % East Sea transport (run)
% end-------------------- get korea strait transport (IPSL-CM5A-LR)


clearvars '*' -except inputyear inputyear1 inputyear2 inputmonth scenname cmip5dir translev



% start-------------------- get korea strait transport (IPSL-CM5A-MR)
testname='IPSL-CM5A-MR';
outputdir=[cmip5dir, 'transport/', scenname, '/Omon/', testname, '/'];
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
	fprintf(fid, '%%Korea  	 \n');
	fclose(fid);
        for monthij=1:length(inputmonth)
            tempmonth=inputmonth(monthij);
            monthstr=num2str(tempmonth, '%02i');
            xData((12*(yearij-1))+monthij) = datenum([num2str(tempmonth,'%02i'),'-01-',num2str(tempyear,'%04i')]);

            varname='vo';
            filedir = strcat(cmip5dir, varname, '/', scenname, '/Omon/', testname, '/'); % % where data files are
            flag_file_in = false;
            list = dir( [ filedir, '/', varname, '*' ]); 
            for kk = 1 : length( list )
                fname_in    = list(kk).name;
                fname_split = strsplit( fname_in, {'_','.'} );
                fyear_str   = strsplit( fname_split{end-1}, '-' );
                fyear_start = str2num( fyear_str{1}(1:4) );
                fyear_end   = str2num( fyear_str{2}(1:4) );
                if( tempyear >= fyear_start && tempyear <= fyear_end &&     ...                 
                        strcmp( fname_split{2}, 'Omon' ) &&         ...
                        strcmp( fname_split{3}, testname ) &&      ...                 
                        strcmp( fname_split{4}, scenname ) )
                    flag_file_in = true;            break;
                end         
            end         
            if( ~flag_file_in )
                fprintf('Source File for %04i does not Exist. Continue...\n',year);   continue;
            end
            vfilename=[filedir, '/', fname_in];
            tind_v=(tempyear-fyear_start)*12+tempmonth;

            v = ncread(vfilename,varname,[27 100 1 tind_v], [1 1 inf 1]);  %% cut horizontal area [x,y,z] (wider than target area)
            lev=ncread(vfilename,'lev');
            diff_lev=diff(lev);
            for thickij=1:length(lev)
                if thickij==1
                    thick_lev(thickij)=lev(thickij)+diff_lev(thickij)/2.0;
                elseif thickij==length(lev)
                    thick_lev(thickij)=diff_lev(thickij-1);
                else
                    thick_lev(thickij)=(diff_lev(thickij-1)+diff_lev(thickij))/2.0;
                end
            end
            vsum=sum(squeeze(v).*thick_lev','omitnan');
            if (exist('lons')==0)
                lons=ncread(vfilename,'lon',[26,100], [1 1]);
                lone=ncread(vfilename,'lon',[28,100], [1 1]);
                lats=ncread(vfilename,'lat',[26,100], [1 1]);
                late=ncread(vfilename,'lat',[28,100], [1 1]);
                xdist=m_lldist([lons,lone], [lats,late]) /2.0 * 1000.0; %% (m)
            end
            korea_tr((yearij-1)*12+monthij)=vsum*xdist/1000000.0;
            disp([testname, num2str(tempyear),'Y', num2str(tempmonth), 'M']);
            fid=fopen(txtname, 'a+');
	    fprintf(fid, '%8.3f \n', korea_tr((yearij-1)*12+monthij));
	    fclose(fid);
        end
    end
% end-------------------- get korea strait transport (IPSL-CM5A-MR)


clearvars '*' -except inputyear inputyear1 inputyear2 inputmonth scenname cmip5dir translev

% start-------------------- get korea strait transport (NorESM1-M)
testname='NorESM1-M';
outputdir=[cmip5dir, 'transport/', scenname, '/Omon/', testname, '/'];
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
	fprintf(fid, '%%Korea  	 \n');
	fclose(fid);
        for monthij=1:length(inputmonth)
            tempmonth=inputmonth(monthij);
            monthstr=num2str(tempmonth, '%02i');
            xData((12*(yearij-1))+monthij) = datenum([num2str(tempmonth,'%02i'),'-01-',num2str(tempyear,'%04i')]);

            varname='vo';
            filedir = strcat(cmip5dir, varname, '/', scenname, '/Omon/', testname, '/'); % % where data files are
            flag_file_in = false;
            list = dir( [ filedir, '/', varname, '*' ]); 
            for kk = 1 : length( list )
                fname_in    = list(kk).name;
                fname_split = strsplit( fname_in, {'_','.'} );
                fyear_str   = strsplit( fname_split{end-1}, '-' );
                fyear_start = str2num( fyear_str{1}(1:4) );
                fyear_end   = str2num( fyear_str{2}(1:4) );
                if( tempyear >= fyear_start && tempyear <= fyear_end &&     ...                 
                        strcmp( fname_split{2}, 'Omon' ) &&         ...
                        strcmp( fname_split{3}, testname ) &&      ...                 
                        strcmp( fname_split{4}, scenname ) )
                    flag_file_in = true;            break;
                end         
            end         
            if( ~flag_file_in )
                fprintf('Source File for %04i does not Exist. Continue...\n',year);   continue;
            end
            vfilename=[filedir, '/', fname_in];
            tind_v=(tempyear-fyear_start)*12+tempmonth;

            v = ncread(vfilename,varname,[151 281 1 tind_v], [2 1 inf 1]);  %% cut horizontal area [x,y,z] (wider than target area)
            lev=ncread(vfilename,'lev');
            diff_lev=diff(lev);
            for thickij=1:length(lev)
                if thickij==1
                    thick_lev(thickij)=lev(thickij)+diff_lev(thickij)/2.0;
                elseif thickij==length(lev)
                    thick_lev(thickij)=diff_lev(thickij-1);
                else
                    thick_lev(thickij)=(diff_lev(thickij-1)+diff_lev(thickij))/2.0;
                end
            end
            thick_lev_2=repmat(thick_lev, [2 1])';
            vsum=sum(sum(squeeze(v).*thick_lev_2','omitnan'),'omitnan');
            if (exist('lons')==0)
                lons=ncread(vfilename,'lon',[150,281], [1 1]);
                lone=ncread(vfilename,'lon',[153,281], [1 1]);
                lats=ncread(vfilename,'lat',[150,281], [1 1]);
                late=ncread(vfilename,'lat',[153,281], [1 1]);
                xdist=m_lldist([lons,lone], [lats,late]) *1.0 /3.0 * 1000.0; %% (m)
            end
            korea_tr((yearij-1)*12+monthij)=vsum*xdist/1000000.0;
            
            disp([testname, num2str(tempyear),'Y', num2str(tempmonth), 'M']);
            fid=fopen(txtname, 'a+');
	    fprintf(fid, '%8.3f \n', korea_tr((yearij-1)*12+monthij));
	    fclose(fid);
        end
    end
% end-------------------- get korea strait transport (NorESM1-M)


clearvars '*' -except inputyear inputyear1 inputyear2 inputmonth scenname cmip5dir translev

% start-------------------- get korea strait transport (MPI-ESM-LR)
testname='MPI-ESM-LR';
outputdir=[cmip5dir, 'transport/', scenname, '/Omon/', testname, '/'];
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
	fprintf(fid, '%%Korea  	 \n');
	fclose(fid);
        for monthij=1:length(inputmonth)
            tempmonth=inputmonth(monthij);
            monthstr=num2str(tempmonth, '%02i');
            xData((12*(yearij-1))+monthij) = datenum([num2str(tempmonth,'%02i'),'-01-',num2str(tempyear,'%04i')]);

            varname='vo';
            filedir = strcat(cmip5dir, varname, '/', scenname, '/Omon/', testname, '/'); % % where data files are
            flag_file_in = false;
            list = dir( [ filedir, '/', varname, '*' ]); 
            for kk = 1 : length( list )
                fname_in    = list(kk).name;
                fname_split = strsplit( fname_in, {'_','.'} );
                fyear_str   = strsplit( fname_split{end-1}, '-' );
                fyear_start = str2num( fyear_str{1}(1:4) );
                fyear_end   = str2num( fyear_str{2}(1:4) );
                if( tempyear >= fyear_start && tempyear <= fyear_end &&     ...                 
                        strcmp( fname_split{2}, 'Omon' ) &&         ...
                        strcmp( fname_split{3}, testname ) &&      ...                 
                        strcmp( fname_split{4}, scenname ) )
                    flag_file_in = true;            break;
                end         
            end         
            if( ~flag_file_in )
                fprintf('Source File for %04i does not Exist. Continue...\n',year);   continue;
            end
            vfilename=[filedir, '/', fname_in];
            tind_v=(tempyear-fyear_start)*12+tempmonth;
            v = ncread(vfilename,varname,[236 114 1 tind_v], [1 1 inf 1]);  %% cut horizontal area [x,y,z] (wider than target area)
            lev=ncread(vfilename,'lev');
            diff_lev=diff(lev);
            for thickij=1:length(lev)
                if thickij==1
                    thick_lev(thickij)=lev(thickij)+diff_lev(thickij)/2.0;
                elseif thickij==length(lev)
                    thick_lev(thickij)=diff_lev(thickij-1);
                else
                    thick_lev(thickij)=(diff_lev(thickij-1)+diff_lev(thickij))/2.0;
                end
            end
            vsum=sum(squeeze(v).*thick_lev','omitnan');
            if (exist('lons')==0)
                lons=ncread(vfilename,'lon',[235,114], [1 1]);
                lone=ncread(vfilename,'lon',[237,114], [1 1]);
                lats=ncread(vfilename,'lat',[235,114], [1 1]);
                late=ncread(vfilename,'lat',[237,114], [1 1]);
                xdist=m_lldist([lons,lone], [lats,late]) /2.0 * 1000.0; %% (m)
            end
            korea_tr((yearij-1)*12+monthij)=vsum*xdist/1000000.0;
            
            disp([testname, num2str(tempyear),'Y', num2str(tempmonth), 'M']);
            fid=fopen(txtname, 'a+');
	    fprintf(fid, '%8.3f \n', korea_tr((yearij-1)*12+monthij));
	    fclose(fid);

        end
    end
% end-------------------- get korea strait transport (MPI-ESM-LR)
