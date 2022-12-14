% %  Created 21-Nov-2022 by Yong-Yub Kim  

clc; clear all; close all;

%% set path
tmp.dropboxpath = '/Volumes/kyy_raid/kimyy/Dropbox';
tmp.fs=filesep;
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);

%% model configuration
dirs.root='/Volumes/kyy_raid/earth.system.predictability/ASSM_EXP';
dirs.yoshi_root='/proj/yoshi/DATA/CESM2_ODA';
dirs.archive=[dirs.root, filesep, 'archive'];
dirs.saveroot='/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/ASSM_EXP';

% config.years=1967:2021;
config.years=1999:2021;

config.months=1:12;
config.scenname='HIST';
config.gridname='f09_g17';
config.obsnames={'projdv7.3', 'en4.2'};
config.ensnames={'ba-10p1', 'ba-10p2', 'ba-10p3', 'ba-10p4', 'ba-10p5', ...
    'ba-20p1', 'ba-20p2', 'ba-20p3', 'ba-20p4', 'ba-20p5'};

config.component='ocn';
config.varnames={'TEMP', 'SALT', 'DIC', 'ALK', 'NO3', 'PO4', 'SiO3'};
% config.varnames={'ALK'};

config.len_t_y = length(config.years);
config.len_t_m = length(config.months);
config.len_t = config.len_t_y * config.len_t_m;
config.len_obs= length(config.obsnames);
config.len_ens= length(config.ensnames);



%% observation configuration (BATS, Excel; .xlsx)

config_obs.staname = 'JMA137E';
% config_obs.sta_lon = [360-66.1690, 360-60.4470]; % degrees west to 0~360 degrees east
% config_obs.sta_lat = [24.7590, 35.6670]; % degrees north
% config_obs.sta_lon = 360-64.1725; % degrees west to 0~360 degrees east
% config_obs.sta_lat = 31.6840; % degrees north

config_obs.avgdepth=[100, 140, 200];

dirs.obsroot = '/Volumes/kyy_raid/kimyy/Observation/JMA/section_137E/obsdata';
% dirs.obssavedir = [dirs.obsroot, filesep, 'mat'];
% mkdir(dirs.obssavedir)

OBS.tmp=0;
for tyi = 1:config.len_t_y
    tmp.yearstr=num2str(config.years(tyi), '%04i');
    for tmi = 1:config.len_t_m
        tmp.monthstr=num2str(config.months(tmi), '%02i');
        dirs.obs_ym=[dirs.obsroot, filesep, tmp.yearstr, tmp.monthstr];
        dirs.obs_WAT=[dirs.obs_ym, filesep, 'WAT-FILE'];
        tmp.list = dir( [ dirs.obs_WAT, filesep, '*' ]);
        for ti=3:length(tmp.list)
            tmp.filename=[dirs.obs_WAT, filesep, tmp.list(ti).name];
%             tmp.data= readtable(tmp.filename, 'FileType', 'text');
            %read file
            tmp.filestr = fileread(tmp.filename);
            %break it into lines
            tmp.filebyline = regexp(tmp.filestr, '\n', 'split');
            %remove empty lines
            tmp.filebyline( cellfun(@isempty,tmp.filebyline) ) = [];
            %split by fields
            tmp.filebyfield = regexp(tmp.filebyline, ',', 'split');
            
            for linei=1:10
                tmp.str=strtrim(tmp.filebyfield{linei}{1});
                if (strcmp(tmp.str(1:4), 'Cast')==1)
                    tmp.line_meta=linei;
                elseif (strcmp(tmp.str(1:4), 'STNN')==1)
                    tmp.line_dat=linei+2; break;
                end
            end
            for wordi=1:length(tmp.filebyfield{tmp.line_meta})
                tmp.str=strtrim(tmp.filebyfield{tmp.line_meta}{wordi});
                if (strcmp(tmp.str, 'Date')==1)
                    tmp.datestr=strtrim(tmp.filebyfield{tmp.line_meta}{wordi+1});
                    tmp.obs_year=str2double(tmp.datestr(1:4));
                    tmp.obs_month=str2double(tmp.datestr(6:7));
                    tmp.obs_day=str2double(tmp.datestr(9:10));
                elseif (strcmp(tmp.str, 'Lat.')==1)
                    tmp.latstr=strtrim(tmp.filebyfield{tmp.line_meta}{wordi+1});
                    tmp.obs_lat=str2double(tmp.latstr(1:2))+str2double(tmp.latstr(4:5))/60.0;
                elseif (strcmp(tmp.str, 'Lon.')==1)
                    tmp.lonstr=strtrim(tmp.filebyfield{tmp.line_meta}{wordi+1});
                    tmp.obs_lon=str2double(tmp.lonstr(1:3))+str2double(tmp.lonstr(5:6))/60.0;
                end
            end
            
%             tmp.line_dat=tmp.filebyfield{1}{1}

            if isfield(OBS, 'data')==0
                OBS.data=tmp.data(2:end,:);
            else
                tmp.tsz=size(tmp.data,1)-1;
                OBS.data(end+1:end+tmp.tsz,:)=tmp.data(2:end,:);
            end
        end
    end
end


config_obs.filename = [dirs.obsroot, filesep, 'BATS/batsftp.bios.edu/BATS/bottle/bats_bottle.xlsx'];


filename='/Volumes/kyy_raid/kimyy/Observation/JMA/section_137E/obsdata/202201/WAT-FILE/ks6090_e4.wat';

fid=fopen(filename);
abc=textscan(fid,'%s');
fclose(fid);
ind_latstr=find(strcmp(string(abc{1}),"Lat.,")==1);
lat=abc{1}{ind_latstr+1};

readtable(filename, 'FileType', 'text');


% %read file
% filestr = fileread(tmp.filename);
% %break it into lines
% filebyline = regexp(filestr, '\n', 'split');
% %remove empty lines
% filebyline( cellfun(@isempty,filebyline) ) = [];
% %split by fields
% filebyfield = regexp(filebyline, ',', 'split');

% strtrim(filebyline{6})

% date -> WAT file list up -> read -> arrangement by time




% filedir = strcat(cmip6dir, variable, '/interp/', testname, filesep, ensname, filesep, 'gn', filesep); % % where data files are
% flag_file_in = false;
% list = dir( [ filedir, filesep, variable, '*' ]); 
% for kk = 1 : length( list )
%     fname_in    = list(kk).name;
%     fname_split = strsplit( fname_in, {'_','.'} );
%     fyear_str   = strsplit( fname_split{end-1}, '-' );
%     fyear= str2num( fyear_str{1}(1:4) );
%     if( tempyear == fyear  &&     ...                 
%             strcmp( fname_split{2}, 'interp' ) &&         ...
%             strcmp( fname_split{3}, testname ) &&      ...                 
%             strcmp( fname_split{4}, scenname ) )
%         flag_file_in = true;            break;
%     end         
% end     