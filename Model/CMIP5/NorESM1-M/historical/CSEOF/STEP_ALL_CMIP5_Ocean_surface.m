clc;clear all;close all;
% % This script is based on MATLAB 2017a
% % Updated 28-Mar-2019 by Y-Y. Kim.

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
% var_names = {'tas', 'psl', 'hur', 'rsds', 'ua', 'va'};  % tas, psl, hur, rsds, ua, va
var_names = {'zos'};  % tas, psl, hur, rsds, ua, va
% var_names = {'thetao'};  % tas, psl, hur, rsds, ua, va

model_name = 'NorESM1-M';
scen_name = 'historical';
inputyear = 1976:2005;
dl = 0.5;
section = [115 164 15 52];  
regress_flag = 1;   % regression switch
tgt_var_name = 'tas';   % target variable name
% if (regress_flag == 1)
%     regress_tgt_mode = 2;
% end

for nvar = 1:length(var_names)
    var_name = var_names{nvar};
    %% ---------case 1. ATM file(ocean mask)-------------%
    %example of temperature (40 layer)
    %-----------STEP1. set the variables----------------%
    %input file path,name
    %path in the damo server
    tfreq = 1;  %% 1 : monthly, 2 : daily

    for nyear = 1:length(inputyear)
      tempyear = inputyear(nyear);
      filepath = [workdir, '/data/', model_name, '_mon/', scen_name, '/', var_name]
      
      filename = [var_name, '_interp_', model_name, '_', ...
          scen_name, '_r1i1p1_', num2str(tempyear,'%04i'), '.nc']
      file = strcat(filepath, '/', filename)

    % % % %  open nc file and set the study area
      if (nyear == 1)
        lon_glo = ncread(file,'lon');
        lat_glo = ncread(file,'lat');
    %    lev = ncread(file,'S');
        xlen_glo=length(lon_glo);
        ylen_glo=length(lat_glo);
    %    zlen=length(lev)
        [indw, inde, inds, indn] = findind_Y(dl, section, lon_glo, lat_glo);

        lon = ncread(file,'lon',indw, inde-indw+1);
        lat = ncread(file,'lat',inds, indn-inds+1);
        xlen = length(lon);
        ylen = length(lat);
      end 
      % first time is 16-JAN-1979, monthly data.
      %time = nc{'TIME1'}(49:408);
      %time = ncread(file, 'TIME1')
      
      predict_file_info=ncinfo(file);
      num_of_var=length(predict_file_info.Variables);
      for i=1:num_of_var
        if (strcmp(predict_file_info.Variables(i).Name,var_name)==1)
            predictor_varind=i;
        end
      end
      num_of_dim=length(predict_file_info.Variables(predictor_varind).Dimensions);
      
      if (tfreq==1)
          tstart=(nyear-1)*12+1;
          tend=nyear*12;
          temptime = ncread(file,'time',1,12);
          if (num_of_dim==3)
            ttemp = ncread(file,var_name,[indw inds 1], [xlen ylen 12]);
          elseif (num_of_dim==4)
            ttemp = ncread(file,var_name,[indw inds 1 1], [xlen ylen 1 12]);
          end
      elseif (tfreq==2)
          tstart=(nyear-1)*365+1;
          tend=nyear*365;
          temptime = ncread(file,'time',1,365);
          ttemp = ncread(file,var_name,[indw inds 1], [xlen ylen 365]);
      end

    %   temptime = ncread(file,'time');
      time(tstart:tend)=temptime;
      tlen=length(time)

      temp(:,:,tstart:tend)=ttemp;
    end

    % j : layer(vertical)
    for j = 1 : 1  
        %setting land mask and save mask file
    %     size(temp2)
        n=length(lon);
        q=length(lat);
        m=length(time);
    %     get salt data for land mask
        saltname = [filepath, '/', var_name, '_interp_', model_name, '_', ...
          scen_name, '_r1i1p1_', num2str(tempyear,'%04i'), '.nc'];
        salt = ncread(saltname, var_name);
        lon_salt =ncread(saltname, 'lon');
        lat_salt =ncread(saltname, 'lat');
        salt_interped=griddata(double(lon_salt),double(lat_salt'), double(squeeze(salt(:,:,1,1))'), double(lon), double(lat'))';
        land_idx = zeros(n*q,1);
        idx = find(isnan(salt_interped)==1) ; %land mask: temp==0
        land_idx(idx) = 1 ;
        land_mask = find(land_idx ==1);
        ocean_mask = find(land_idx ==0);
        ocean_mask_name=strcat(workdir, '/data/', model_name, '_mon/ocean_mask_layer_', num2str(j,'%02i'), '.txt');
        save(ocean_mask_name,'ocean_mask','-ascii');

        %extract temp data(NaN data remove)
        data =zeros(n*q- length(land_mask),m);
        for i = 1:m 
            temp2=squeeze(temp(:,:,i));
        %temp2(i,j,:,:)=permute(temp(:,:,j,i), [4 3 2 1]); %temp : lon, lat , depth, time
            temp2(land_mask) = [];
            data(:,i) = temp2(:)' ; % [space time]
    %	size(temp2)
        end

        sizedata=size(data); 
        %% for fast calculation speed, if [M*N], M must be higher than N.
        if (sizedata(1) >= sizedata(2)) 
            data2 = data;
        else
            data2 = data';
        end
        data_len = n*q-length(land_mask)
        size(data2)
        slen = size(data2,1);
        %save data
        tt = strcat(workdir, '/data/', model_name, '_mon/', scen_name, '/', var_name, '/', model_name, '_', var_name, '_', num2str(j,'%02i'), '.data');
        save(tt,'data2','-ASCII');
        disp(['dimension of sampling stations : ',num2str(slen)]);
        disp(['number of sampling points at each station : ',num2str(tlen)]);
    end



    %% STEP 2_1 Makescript_eigen

    % make cseof directory for each variable, model
%     var_name='hur'
    cd([workdir, '/script/']);
    cseof_output_dir = [workdir, '/cseofs/', model_name, '/', scen_name, '/', var_name]
    mkdir(cseof_output_dir);

    name= cell(28);
    %----------STEP2. make script file-----------------------------------%
    name{1} = '#!/bin/csh';     %use csh shell
    % for i = 1:length(nn);
    name{2} = char('');
    name{3} = char(['gfortran -fno-backtrace -o eigen ', workdir, '/programs/eigen/eigenx.f']); %compile the eigen file
    name{4} = char('');
    name{5} = char('cat >! eof.com <<ENDc'  );
%     name1 = char(name1, strcat('../data/',nn(i,:),'.data')  );         %load variable file
%         name{6} = char(tt);  %load variable file. if it is longer than 50 character, you must correct filename variable size in eigenx.f or eigen.f
    ttt = ['./../data/', model_name, '_', var_name, '_', num2str(j,'%02i'), '.data'];
    system(['ln -sf ', tt, ' ', ttt]);  %% eigen fortran codes cannot recognize filenames longer than 50 characters.
    name{6} = char(ttt);
    name{7} = char(['((', num2str(tlen), 'e16.7))'] );  %ascii file read, (time_num)e16.7
    name{8} = char([num2str(slen), ' 1']); %space number
    name{9} = char(num2str(tlen));  %time number
    name{10} = char('1            '    );  %space*time matrix or time*space matrix in matlab(1 or 2)
    name{11} = char('0            ' );                            
    name{12} = char('1          '  );
    name{13} = char('0              ' ); %% area adjustment (0: No)
    name{14} = char('15. 37.' );  %% starting latitude and increment for area adjustment
    name{15} = char('99.99');  %want to modes to explain 95% (percent variance)
    name{16} = char('1.' );  %% EOF scaling factor
    eofmodenum= 10;
    maxmodenum = num2str(eofmodenum);
    name{17} = char(maxmodenum);  %% number of EOFs to be printed                              
    name{18} = char('0');  %% PC normalization
%     name{} = char(name{}, strcat('../cseofs/eof_',nn(i,:),'.dat')   );  %save eof LV file
%         name{19} = char(strcat(cseof_output_dir, '/eof_', model_name, '_', var_name, '.dat')   );  %save eof LV file    
    eof_tt = strcat(cseof_output_dir, '/eof_', model_name, '_', var_name, '.dat');
    eof_ttt = strcat('./../data/eof_', model_name, '_', var_name, '.dat');
    name{19} = char(eof_ttt);  %save eof LV file    
    name{20} = char('DIR');                                       %DIR : binary
%     name{} = char(name{}, strcat('../cseofs/pct_',nn(i,:),'.dat')   );  %save eof pct file
    pct_tt = strcat(cseof_output_dir, '/pct_', model_name, '_', var_name, '.dat');
    pct_ttt = strcat('./../cseofs/pct_', model_name, '_', var_name, '.dat');    
    name{21} = char(pct_ttt);  %save eof pct file
    name{22} = char('(5e13.5)' );     %% ASCII format "in parenthesis"                               %save option: 6e13.5 acsii
    name{23} = strtrim(char('ENDc'));
    name{24} = char(['./eigen < ./eof.com'] );
%     name{} = char(name{}, [workdir, '/script/eigen < ', workdir, '/script/eof.com'] );
    name{25} = char(['mv -f ', workdir, '/script/inform.d ', cseof_output_dir, '/inf_', 'NorESM1-M_', var_name,'.d']); %save information file
    name{26} = char(['mv -f ', workdir, '/script/avg.d ', cseof_output_dir, '/avg_', 'NorESM1-M_', var_name,'.d']);  %save average file
    name{27} = char('rm -f eigen eof.com' );
    name{28} = char('' );
    %     name = char(name,strtrim(name1));
    % end

    fid = fopen([workdir, '/script/eigen.c'], 'w+')
    for nline=1:length(name)
        fprintf(fid, '%s\n', name{nline});
    end
    fclose(fid);
    system(['cd ', workdir, '/script'])
    system(['csh -xv eigen.c > ', cseof_output_dir, '/eigen_log.log']);
        
    %% STEP 2_1 Makescript_cseof
    name = cell(28);
    name{1} = '#!/bin/csh';     %use csh shell
    % for i = 1:length(nn);
    name{2} = char('');
    name{3} = char( ['gfortran -fno-backtrace -o cseof ', workdir, '/programs/eigen/cseof.f']); %compile the eigen file
    name{4} = char( '' );
    name{5} = char( 'cat >! cseof.com <<ENDc'  );
    name{6} = char('0');  %% job number (only 0)
    name{7} = char( pct_ttt  );     %load variable file
    name{8} = char( '(5e13.5)' );
    name{9} = char( [maxmodenum, ' 1']  ); %max.mode number(same to eigen.c)
    name{10} = char( num2str(tlen)); %time number
    name{11} = char(  '1'    ); %% time index(=1) or space index(=2) first
    LVnumber = length(time)/length(inputyear);
    name{12} = char(  num2str(LVnumber) );   %want period(LV number)
    spectralpoints = floor(LVnumber/2);
    name{13} = char(  num2str(spectralpoints)  );   %want period/2 (number of spectral points)
    name{14} = char( '1              ' );  %% interval subdivisions for integrations
    name{15} = char( '0' );  %% cycle for detrending
    name{16} = char( num2str(tlen));    %time number (size of cov matrix)
    name{17} = char( '99.99' );   %want explain %
    cseofmodenum=10;
    name{18} = char( num2str(cseofmodenum) );    %want max.mode number
    name{19} = char( '1.');   %% eof scale factor (Do not change)
    name{20} = char( '2');  %% 1: rc ts    2 : cov
    name{21} = char( 'ENDc' );
    name{22} = char( [workdir, '/script/cseof < ', workdir, '/script/cseof.com'] );
    cinf_tt = strcat(cseof_output_dir, '/cinf_', model_name, '_', var_name, '.d');
    cinf_ttt = strcat('./../cseofs/cinf_', model_name, '_', var_name, '.d');
    blo_tt = strcat(cseof_output_dir, '/blo_', model_name, '_', var_name, '.d');
    blo_ttt = strcat('./../cseofs/blo_', model_name, '_', var_name, '.d');
    cpct_tt = strcat(cseof_output_dir, '/cpct_', model_name, '_', var_name, '.d');
    cpct_ttt = strcat('./../cseofs/cpct_', model_name, '_', var_name, '.d');
    ceig_tt = strcat(cseof_output_dir, '/ceig_', model_name, '_', var_name, '.d');
    ceig_ttt = strcat('./../cseofs/ceig_', model_name, '_', var_name, '.d');
    name{23} = char( ['mv -f inform.d ', cinf_ttt]);
    name{24} = char(  ['mv -f Bloch.d ', blo_ttt] );
    name{25} = char(  ['mv -f pcts.d ', cpct_ttt] );
    name{26} = char(  ['mv -f eigen.d ', ceig_ttt] );
    name{27} = char( 'rm -f cseof cseof.com' );
    name{28} = char( '' );
    % end

    fid = fopen([workdir, '/script/cseof.c'], 'w+')
    for nline=1:size(name,1)
        fprintf(fid, '%s\n', name{nline});
    end
    fclose(fid);
    system(['csh -xv cseof.c > ', cseof_output_dir, '/cseof_log.log']);
    
    
    
    
    %% STEP 2_5 Makescript_recast
    %------case1. one layer-----------------%
    name = cell(23);
    name{1} = '#!/bin/csh';     %use csh shell

    % for i = 2:2;
    name{2} = char('');
    name{3} = char('gfortran -fno-backtrace -o recast ../programs/util/recastx.f');
    name{4} = char( '' );
    name{5} = char( 'cat >! recast.com <<ENDc'  );
    name{6} = char( eof_ttt  ); %predictor eof file
    name{7} = char( 'DIR' );
    name{8} = char( 'nofile'  );
    name{9} = char( num2str(cseofmodenum));                     %max. eof mode number
    name{10} = char(  num2str(slen) );   %space number
    name{11} = char(  '1.      ' );                              
    name{12} = char(  [num2str(slen), ' 1']  );     %space structure
    % name{13} = char( strcat('../regress/blo_reg_', 'NorESM1-M_', varname,'.d') ); %regressed blo file
    name{13} = char( blo_ttt ); %blo file
    name{14} = char( 'SEQ' );
    name{15} = char( '0');
    sizeLV = LVnumber * cseofmodenum;
    name{16} = char( num2str(sizeLV) );  %LV number * cseof mode number (365*10)
    name{17} = char( '1' );
    LV_tt = strcat(cseof_output_dir, '/LV_', model_name, '_', var_name, '.dat');
    LV_ttt = strcat('./../cseofs/LV_', model_name, '_', var_name, '.dat');
    name{18} = char( LV_ttt);     %LV file name
    name{19} = char( 'DIR');
    name{20} = char( 'ENDc' );
    name{21} = char( [workdir, '/script/recast < ', workdir, '/script/recast.com'] );
    name{22} = char( 'rm -f recast recast.com' );
    name{23} = char( '' );
    % end

    fid = fopen([workdir, '/script/recast.c'], 'w+')
    for nline=1:size(name,1)
        fprintf(fid, '%s\n', name{nline});
    end
    fclose(fid);
    system(['csh -xv recast.c > ', cseof_output_dir, '/recast_log.log']);
    
    system(['mv -f ', eof_ttt, ' ', eof_tt]); 
    system(['mv -f ', pct_ttt, ' ', pct_tt]);  
    system(['mv -f ', cinf_ttt, ' ', cinf_tt]); 
    system(['mv -f ', blo_ttt, ' ', blo_tt]);  
    system(['mv -f ', cpct_ttt, ' ', cpct_tt]); 
    system(['mv -f ', ceig_ttt, ' ', ceig_tt]);  
    system(['mv -f ', LV_ttt, ' ', LV_tt]);  

    system(['rm -f ', ttt]);
    
    cd(presentdir);
    
    
    %% STEP 3 LV_layer_merge
    filepath = [workdir, '/cseofs/', model_name, '_mon/', scen_name, '/', var_name]
    
    filename = [var_name, '_interp_', model_name, '_', scen_name, ...
        '_', num2str(tempyear,'%04i'), '.nc']

    file = strcat(filepath, '/', filename)

    %----------STEP1. load lon,lat information-------------%
    %load atm lon,lat information
    lon_atm = lon;
    lat_atm = lat;
    lon_uv = lon;
    lat_uv = lat;
    lon2 = lon;
    lat2 = lat;

    %% ocean component

    %thetao
    T= LVnumber;
    mode = cseofmodenum;
    LV_increment1 = zeros(1,length(lon2),length(lat2),T,mode);
    for layer = 1 : 1
        %set filepath , variables
        filepath = [workdir, '/cseofs/', model_name, '/', scen_name, '/', var_name]
        LV_name = strcat(filepath,'/', 'LV_NorESM1-M_', var_name, '.dat');
        %load mask data
        ocean_mask = importdata(ocean_mask_name);
        time = size(data,2) ; % time(month)
        T = LVnumber  ;  %period
    %     mode = str2num(mode) ; %number of LV mode
        lon = lon2; %lon
        lat=  lat2; %lat
        mm = length(lon);
        nn = length(lat);

        %%%%Loading vector
        %load LV data and reshape
        raw_LV = fopen(LV_name);
        raw_LV = fread(raw_LV,'float');
        space = length(raw_LV)/mode/T;
        raw_LV = reshape(raw_LV,space,T*mode);
        LV_grid = zeros(nn*mm,T*mode) ;
        LV_grid(ocean_mask,:) = raw_LV(:,:);
        LV = reshape(LV_grid,mm,nn,T,mode);
        LV_increment1(layer,:,:,:,:) = LV;
    end
    save([cseof_output_dir, '/lv_layer_merge.mat'], '-v7.3');
    
% % %     STEP4, PLOT
    rehash toolboxcache
    figure_output_dir = [workdir, '/figure/', model_name, '/', scen_name, '/', var_name]
    mkdir([figure_output_dir]);
% % %     %load loading vector file
% % %     load([cseof_output_dir, '/lv_layer_merge.mat']);

    mode = 3;
    %load cseof pc time series
    pct_data = importdata(cpct_tt)';
    pct = pct_data(:)';
    pct = reshape(pct,tlen,length(pct)/tlen);
    %pct(:,2) = -pct(:,2);
    %plot pc time series

    for nyear = 1:length(inputyear)
        tempyear = inputyear(nyear);
        for month=1:12
            xData((12*(nyear-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
        end
    end
    
    clear i
    for i = 1:mode
        plot(xData, pct(:,i));
        tt = [var_name, '-',num2str(i, '%02i'),'mode'];
        title(tt,'fontsize',18,'fontweight','bold');  
        datetick('x', 'yy', 'keepticks')
        xlabel('Time(year)','FontSize',20) ;
        set(gca,'FontSize',14);
        saveas(gcf, strcat([figure_output_dir, '/', tt, '.png']),'png');
        close all
        i
    end

    %LV = LV(layer,lon,lat,period,mode)
    LV_increment1(LV_increment1==0) = NaN;
    %LV_thetao(:,:,:,:,2) = -LV_thetao(:,:,:,:,2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%% surface plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     hh = 1;
    %plot - surface
    [lon2, lat2] = meshgrid(lon,lat);
    %[lat2, lon2] = meshgrid(lat,lon);
    
    mkdir([figure_output_dir, '/surface']);
    %sst

    for i = 1 : mode
        for period = 1:T
        %dat3 = squeeze(LV_increment1(hh,:,:,period,i));
        dat3 = squeeze(LV_increment1(hh,:,:,period,i));
        dat3(dat3==1e20) = NaN;
        m_proj('mercator','lon',[115 164],'lat',[15 52]);
        m_gshhs_l('color','k');
        m_gshhs_l('patch',[.8,.8,.8]);
        m_grid('box','fancy','tickdir','in','linewidth',1);
        hold on;
        %[C,h] = m_contour(lon2,lat2,dat3,'LineWidth',2);
        m_pcolor(lon2,lat2,dat3');
        %clabel(C,h,'FontSize',14,'Color','k','Rotation',0,'fontweight','bold');
        shading interp;
        colorbar;
%         colormap jet;
        bwr_map
        colormap(bwrmap);
        xlabel(['Longitude (^o E)'],'fontsize',18,'fontweight','bold','fontname','times new roman');
        ylabel(['Latitude(^o N)'],'fontsize',18,'fontweight','bold','fontname','times new roman');
        tt = [var_name, '-',num2str(i,'%02i'),'mode','(period-',num2str(period,'%02i'),')'];
        title(tt,'fontsize',18,'fontweight','bold');   %��ٲٱ�
%         caxis([-4 4]);
        clim_val=max(max(abs(dat3')));
        caxis([-clim_val clim_val]);
        axis tight;
        saveas(gcf,strcat([figure_output_dir, '/surface/',tt,'.png']),'png');
        close all;
        [i period]
        end
    end
    
% % % %% STEP 2_3 Makescript_regress
    if (regress_flag == 1 && strcmp(tgt_var_name,var_name) ~= 1)
        
        cd([workdir, '/script/']);

        tgt_cseof_output_dir = [workdir, '/cseofs/', model_name, '/', scen_name, '/', tgt_var_name];
        tgt_cpct_tt = strcat(tgt_cseof_output_dir, '/cpct_', model_name, '_', tgt_var_name, '.d');
        tgt_cpct_ttt = strcat('./../cseofs/cpct_', model_name, '_', tgt_var_name, '.d');
        tgt_cinf_tt = strcat(tgt_cseof_output_dir, '/cinf_', model_name, '_', tgt_var_name, '.d');
        tgt_cinf_ttt = strcat('./../cseofs/cinf_', model_name, '_', tgt_var_name, '.d');
        
        system(['ln -sf ', cinf_tt, ' ', cinf_ttt]); 
        system(['ln -sf ', cpct_tt, ' ', cpct_ttt]); 
        system(['ln -sf ', tgt_cinf_tt, ' ', tgt_cinf_ttt]); 
        system(['ln -sf ', tgt_cpct_tt, ' ', tgt_cpct_ttt]); 
        
        for regress_tgt_mode=1:cseofmodenum
            %----This is the code to make regress.c script file for same grid but
            %differnet variables

            %----------STEP2. make script file-----------------------------------%
            
            name = cell(30);
            name{1} = '#!/bin/csh'; %use csh shell
            name{2} =  char('');
            name{3} =  char('gfortran -o regress ../programs/util/regress_new_mean.f');
            name{4} =  char( '' );
            name{5} =  char( ['cat >! regress', num2str(regress_tgt_mode, '%02i'), '.com << END']  );
            name{6} =  char([num2str(tlen), ' 10    ']);   %Want regress PC time series (time,predicted mode)
            name{7} =  char( tgt_cpct_ttt );        %target PC file name
            name{8} =  char( '(6e13.5)' );                     %format
            regress_skip_line=(regress_tgt_mode-1)*tlen/6;
            name{9} =  char( num2str(regress_skip_line)  );                 %skip amount   
            name{10} =  char( ['1 ', num2str(tlen)]);                    %data interval for target time series
            name{11} =  char(  cpct_ttt    ); %predictor time series
            name{12} =  char(  '(6e13.5)    ' );                     %format
%             name{13} =  char( num2str(regress_skip_line)  );                 %skip amount   
            name{13} =  char( '0'  );                     %skip amount
            name{14} =  char( ['1 ', num2str(tlen)] );    %data interval for predictor time series
            name{15} =  char( ['1 ', num2str(tlen)] );    %regression interval
            name{16} =  char( 'regress.d');
            name{17} =  char( '1         ' );  %estimated output (0: No,  1: Yes)
            name{18} =  char( '2        ' );   %confidence interval(0:no,2:90%,3:95%,4:99%)
            name{19} =  char( '2          ');  % scaling option (just use 2)
            name{20} =  char( tgt_cinf_ttt);    %target cinf file
            if (regress_tgt_mode==1)
                name{21} =  char( '(41x,e15.7)' );                             %format
                num_slash=['/'];
            else
                num_slash=[num_slash,'/'];
                format_tgt = ['(', num_slash, ',41x,e15.7)'];
                name{21} =  char( format_tgt );                             %format
            end
            name{22} =  char( cinf_ttt );   %predictor cinf file
            name{23} =  char( '(41x,e15.7) ' );                   %predictor format. put first line's format (because program reads variables of all modes)
            name{24} =  char( 'END' );
            name{25} =  char( ['./regress < regress', num2str(regress_tgt_mode, '%02i'),'.com'] );
            name{26} =  char( 'mv regress.d regress.s');
            name{27} =  char( '' ); 
            %save script
            reg_tt = strcat(cseof_output_dir, '/reg_', model_name, '_', var_name, '.d');
            reg_ttt = strcat('./../cseofs/reg_', model_name, '_', var_name, '.d');
            if (regress_tgt_mode==1)
                name{28} = char( ['mv -f regress.s ', reg_ttt] );
            else
                name{28} = char( ['cat regress.s >> ', reg_ttt] );
            end
            name{29} = char('rm -f regress regress*.com regress.d regress.s' );
            name{30} = char('rm -f Y_err.d Y_est.d' );

            fid = fopen([workdir, '/script/regress.c'], 'w+')
            for nline=1:size(name,1)
                fprintf(fid, '%s\n', name{nline});
            end
            fclose(fid);

            system(['csh -xv regress.c > ', cseof_output_dir, '/regress_log.log']);
             
        end % regress_target_mode
        
        % % % %% STEP 2_4 Makescript_regress_combine
        
        system(['ln -sf ', blo_tt, ' ', blo_ttt]); 
            
        %----------STEP2. make script file-----------------------------------%
        name = cell(24);

        name{1} = '#!/bin/csh';     %use csh shell
        name{2} = char('');
        name{3} = char('gfortran -o combin ../programs/util/combin_new.f');
        name{4} = char( '' );
        name{5} = char( 'cat >! combin.com << ENDc '  );
        name{6} = char( blo_ttt  );     %predictor blo file
        name{7} = char( 'SEQ' );
        name{8} = char( cinf_ttt );    %predictor inf file
        name{9} = char( ['(22x,e15.7,////,',num2str(cseofmodenum),'(41x,e15.7,/))']  );        %format, 10 should be changed to regress mode number
        name{10} = char( num2str(cseofmodenum));   %regressed mode number
        name{11} = char(  num2str(LVnumber)    );                       %LV number
        name{12} = char(  [num2str(cseofmodenum),' 1          ']  );   %max.EOF mode number of predictor :20(?)
        name{13} = char( num2str(cseofmodenum) );  %regressed mode number. total number of target PCs
        name{14} = char( '0' );                  
        reg_blo_tt = strcat(cseof_output_dir, '/reg_blo_', model_name, '_', var_name, '.d');
        reg_blo_ttt = strcat('./../cseofs/reg_blo_', model_name, '_', var_name, '.d');
        name{15} = char( reg_blo_ttt);  %predictor,regressed blo file
        name{16} = char( 'SEQ' );
        name{17} = char( '2' );
        name{18} = char( reg_ttt);      %target-predictor regress file name
        name{19} = char( ['(///,', num2str(cseofmodenum), '(7X,E13.5,/),/)']);   %10 is regress mode number
        name{20} = char( tgt_cinf_ttt);                     %target cinf file name
        name{21} = char( ['(22x,e15.7,////,', num2str(cseofmodenum),'(41x,e15.7,/))']);  %10 is target mode number
        name{22} = char( 'ENDc' );
        name{23} = char( './combin < combin.com  ' );
        name{24} = char( 'rm -f combin combin.com  ' );

        fid = fopen([workdir, '/script/combin.c'], 'w+')
        for nline=1:size(name,1)
            fprintf(fid, '%s\n', name{nline});
        end
        fclose(fid);
        system(['csh -xv combin.c > ', cseof_output_dir, '/regress_combin_log.log']);


        
        
        % % % %% STEP 2_5 Makescript_regress_recast
        
        system(['ln -sf ', eof_tt, ' ', eof_ttt]); 

        %------case1. one layer-----------------%
        name = cell(23);
        name{1} = '#!/bin/csh';     %use csh shell

        % for i = 2:2;
        name{2} = char('');
        name{3} = char('gfortran -fno-backtrace -o recast ../programs/util/recastx.f');
        name{4} = char( '' );
        name{5} = char( 'cat >! recast.com <<ENDc'  );
        name{6} = char( eof_ttt  ); %predictor eof file
        name{7} = char( 'DIR' );
        name{8} = char( 'nofile'  );
        name{9} = char( num2str(cseofmodenum));                     %max. eof mode number
        name{10} = char(  num2str(slen) );   %space number
        name{11} = char(  '1.      ' );                              
        name{12} = char(  [num2str(slen), ' 1']  );     %space structure
        name{13} = char( reg_blo_ttt ); %regressed blo file
        name{14} = char( 'SEQ' );
        name{15} = char( '0');
        sizeLV = LVnumber * cseofmodenum;
        name{16} = char( num2str(sizeLV) );  %LV number * cseof mode number (365*10)
        name{17} = char( '1' );
        reg_LV_tt = strcat(cseof_output_dir, '/reg_LV_', model_name, '_', var_name, '.d');
        reg_LV_ttt = strcat('./../cseofs/reg_LV_', model_name, '_', var_name, '.d');
        name{18} = char( reg_LV_ttt);     %LV file name
        name{19} = char( 'DIR');
        name{20} = char( 'ENDc' );
        name{21} = char( [workdir, '/script/recast < ', workdir, '/script/recast.com'] );
        name{22} = char( 'rm -f recast recast.com' );
        name{23} = char( '' );
        % end

        fid = fopen([workdir, '/script/recast.c'], 'w+')
        for nline=1:size(name,1)
            fprintf(fid, '%s\n', name{nline});
        end
        fclose(fid);
        
        system(['csh -xv recast.c > ', cseof_output_dir, '/regress_recast_log.log']);
        
        system(['rm -f ', cinf_ttt]); 
        system(['rm -f ', cpct_ttt]); 
        system(['rm -f ', tgt_cpct_ttt]); 
        system(['rm -f ', tgt_cinf_ttt]); 
        system(['rm -f ', blo_ttt]); 
        system(['mv -f ', reg_ttt, ' ', reg_tt]); 
        system(['mv -f ', reg_blo_ttt, ' ', reg_blo_tt]); 
        system(['mv -f ', reg_LV_ttt, ' ', reg_LV_tt]); 

        
        cd(presentdir);
        
        
% % %         combin CPCT from regression coefficient and predictor's PCT
        fid=fopen([cseof_output_dir,'/reg_', model_name, '_', var_name, '.d']);
        for regmode=1:cseofmodenum
            uselessline=fgetl(fid);
            uselessline=fgetl(fid);
            uselessline=fgetl(fid);
            for predictormode=1:cseofmodenum
                abc=fgetl(fid);
                reg_coef(regmode,predictormode)=str2num(abc(9:20));
            end
            uselessline=fgetl(fid);
            abc=fgetl(fid);
            reg_r_sq(regmode,1)=str2num(abc(13:24));
            reg_r_sq(regmode,2)=str2num(abc(26:37));
        end
        fclose(fid);
        
        pct_data = importdata(cpct_tt)';   %% need to change
        pct = pct_data(:)';
        pct = reshape(pct,tlen,length(pct)/tlen);
        
        reg_cpct=zeros(tlen,regmode);
        for regmode=1:cseofmodenum
            for predictormode=1:cseofmodenum
                reg_cpct(:,regmode)=reg_cpct(:,regmode) + pct(:,predictormode) .* reg_coef(regmode,predictormode);
            end
        end
%         plot(xData,reg_cpct(:,1))
        reg_cpct_1line = reshape(reg_cpct,[tlen*regmode 1]);
        
        reg_cpct_tt=[cseof_output_dir,'/reg_cpct_', model_name, '_', var_name, '.d'];
        fid=fopen(reg_cpct_tt, 'w');
        fprintf(fid, '%13.5e%13.5e%13.5e%13.5e%13.5e%13.5e\n', reg_cpct_1line);
        fclose(fid);
        
        %% STEP 3 LV_layer_merge
        filepath = [workdir, '/cseofs/', model_name, '_mon/', scen_name, '/', var_name]
        filename = [var_name, '_interp_', model_name, '_', scen_name, ...
            '_', num2str(tempyear,'%04i'), '.nc']

        file = strcat(filepath, '/', filename)

        %----------STEP1. load lon,lat information-------------%
        %load atm lon,lat information
        lon_atm = lon;
        lat_atm = lat;
        lon_uv = lon;
        lat_uv = lat;
        lon2 = lon;
        lat2 = lat;

        %% ocean component

        %thetao
        T= LVnumber;
        mode = cseofmodenum;
        LV_increment1 = zeros(1,length(lon2),length(lat2),T,mode);
        for layer = 1 : 1
            %set filepath , variables
            filepath = [workdir, '/cseofs/', model_name, '/', scen_name, '/', var_name]
            LV_name = reg_LV_tt;
            %load mask data
            ocean_mask = importdata(ocean_mask_name);
            time = size(data,2) ; % time(month)
            T = LVnumber  ;  %period
        %     mode = str2num(mode) ; %number of LV mode
            lon = lon2; %lon
            lat=  lat2; %lat
            mm = length(lon);
            nn = length(lat);

            %%%%Loading vector
            %load LV data and reshape
            raw_LV = fopen(LV_name);
            raw_LV = fread(raw_LV,'float');
            space = length(raw_LV)/mode/T;
            raw_LV = reshape(raw_LV,space,T*mode);
            LV_grid = zeros(nn*mm,T*mode) ;
            LV_grid(ocean_mask,:) = raw_LV(:,:);
            LV = reshape(LV_grid,mm,nn,T,mode);
            LV_increment1(layer,:,:,:,:) = LV;
        end
        save([cseof_output_dir, '/reg_lv_layer_merge.mat'], '-v7.3');

    % % %     STEP4, PLOT
        rehash toolboxcache
        figure_output_dir = [workdir, '/figure/', model_name, '/', scen_name, '/', var_name]
        mkdir([figure_output_dir]);
    % % %     %load loading vector file
    %     load([cseof_output_dir, '/lv_layer_merge.mat']);

        mode = 3;
        %load cseof pc time series
        pct_data = importdata(reg_cpct_tt)';   %% need to change
        pct = pct_data(:)';
        pct = reshape(pct,tlen,length(pct)/tlen);
        %pct(:,2) = -pct(:,2);
        %plot pc time series

        for nyear = 1:length(inputyear)
            tempyear = inputyear(nyear);
            for month=1:12
                xData((12*(nyear-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
            end
        end

        for i = 1:mode
            plot(xData, pct(:,i));
            tt = ['reg-', var_name, '-',num2str(i,'%02i'),'mode'];
            title(tt,'fontsize',18,'fontweight','bold');  
            datetick('x', 'yy', 'keepticks')
            xlabel('Time(year)','FontSize',20) ;
            set(gca,'FontSize',14);
            saveas(gcf, strcat([figure_output_dir, '/', tt, '.png']),'png');
            close all
            i
        end

        %LV = LV(layer,lon,lat,period,mode)
        LV_increment1(LV_increment1==0) = NaN;
        %LV_thetao(:,:,:,:,2) = -LV_thetao(:,:,:,:,2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%% surface plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         hh = 1;
        %plot - surface
        [lon2, lat2] = meshgrid(lon,lat);
        %[lat2, lon2] = meshgrid(lat,lon);

        mkdir([figure_output_dir, '/surface']);
        %sst
        for i = 1 : mode
            for period = 1:T
            %dat3 = squeeze(LV_increment1(hh,:,:,period,i));
            dat3 = squeeze(LV_increment1(hh,:,:,period,i));
            dat3(dat3==1e20) = NaN;
            m_proj('mercator','lon',[115 164],'lat',[15 52]);
            m_gshhs_l('color','k');
            m_gshhs_l('patch',[.8,.8,.8]);
            m_grid('box','fancy','tickdir','in','linewidth',1);
            hold on;
            %[C,h] = m_contour(lon2,lat2,dat3,'LineWidth',2);
            m_pcolor(lon2,lat2,dat3');
            %clabel(C,h,'FontSize',14,'Color','k','Rotation',0,'fontweight','bold');
            shading interp;
            colorbar;
%             colormap jet;
            bwr_map
            colormap(bwrmap);
            xlabel(['Longitude (^o E)'],'fontsize',18,'fontweight','bold','fontname','times new roman');
            ylabel(['Latitude(^o N)'],'fontsize',18,'fontweight','bold','fontname','times new roman');
            tt = ['reg-', var_name, '-',num2str(i,'%02i'),'mode','(period-',num2str(period,'%02i'),')'];
            title(tt,'fontsize',18,'fontweight','bold');   %��ٲٱ�
%             caxis([-4 4]);
            clim_val=max(max(abs(dat3')));
            caxis([-clim_val clim_val]);
            axis tight;
            saveas(gcf,strcat([figure_output_dir, '/surface/',tt,'.png']),'png');
            close all;
            [i period]
            end
        end
    end % if regress flag 
    
%     nc_tt = strcat(cseof_output_dir, '/cseofs_', model_name, '_', var_name, '.nc');
% 
%     ncid = netcdf.create(nc_tt, 'CLOBBER');
%     xi_rho_dimid = netcdf.defDim(ncid, 'xi_rho', Lp);
%     eta_rho_dimid = netcdf.defDim(ncid, 'eta_rho', Mp);
%     time_dimid = netcdf.defDim(ncid, 'time', 0);
%     eofmode_dimid = netcdf.defDim(ncid, 'eof_mode', 0);
%     cseofmode_dimid = netcdf.defDim(ncid, 'cseof_mode', 0);
%     
%     time_varid = netcdf.defVar(ncid, 'time', 'double', [time_dimid]);
%     netcdf.putAtt(ncid, time_varid, 'long_name', 'day since 0000-01-00 00:00:00');
%     netcdf.putAtt(ncid, time_varid, 'units', 'day');
%     
%     eof_LV_varid = netcdf.defVar(ncid, '');
    
end  % for var_name