clc;clear all;close all;
% % This script is based on MATLAB 2017a
% % Updated 28-Mar-2019 by Y-Y. Kim.
% % Updated 25-Sep-2024 by Y-Y. Kim.

warning off

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
elseif (strcmp(system_name,'MACI64'))
    dropboxpath='/Volumes/kyy_raid/kimyy/Dropbox';
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
end

presentdir = pwd;

workdir='/Users/kimyy/CSEOF'; % %% eigen fortran codes cannot recognize filenames longer than 50 characters. must be short.
var_names = {'sst'};  % tas, psl, hur, rsds, ua, va ...

model_name = 'ERSST';
scen_name = 'historical';
inputyear = 1981:2020;

regress_flag = 0;   % regression switch
tgt_var_name = 'tas';   % target variable name
% if (regress_flag == 1)
%     regress_tgt_mode = 2;
% end




for nvar = 1:length(var_names)
    var_name = var_names{nvar};
    %% ---------case 1. ERSST file (ocean mask)-------------%
    %example of ERSST
    %-----------STEP1. set the variables----------------%
    %input file path,name
    %path in the damo server
    tfreq = 1;  %% 1 : monthly, 2 : daily

    for nyear = 1:length(inputyear)
        tempyear = inputyear(nyear);
        for nmon = 1:12
          filepath = [workdir, '/NWP'];
          filename = ['NWP_ersst_reg_cesm2.v5.', num2str(tempyear,'%04i'), num2str(nmon, '%02i'), '.nc'];
          file = strcat(filepath, '/', filename);
        
        %%  open nc file
          if (nyear == 1)
            lon = ncread(file, 'lon');
            lat = ncread(file, 'lat');
    
            xlen = length(lon);
            ylen = length(lat);
          end 
          
          predict_file_info=ncinfo(file);
          num_of_var=length(predict_file_info.Variables);
          for i=1:num_of_var
            if (strcmp(predict_file_info.Variables(i).Name,var_name)==1)
                predictor_varind=i;
            end
          end
          num_of_dim=length(predict_file_info.Variables(predictor_varind).Dimensions);
          
          if (tfreq==1)
            temptime = ncread(file,'time');
            ttemp = ncread(file,var_name);
          elseif (tfreq==2)
              tstart=(nyear-1)*365+1;
              tend=nyear*365;
              temptime = ncread(file,'time',1,365);
              ttemp = ncread(file,var_name,[indw inds 1], [xlen ylen 365]);
          end
    
          time((nyear-1)*12+nmon)=temptime;
          temp(:,:,(nyear-1)*12+nmon)=ttemp;
        end
    end
    tlen=length(time);

    % j : layer(vertical)
    for j = 1 : 1  
        %% setting land&ocean mask and save mask file for recast
    %     size(temp2)
        n=length(lon);
        q=length(lat);
        m=length(time);
        temp_mask=squeeze(temp(:,:,1));
        land_idx = zeros(n*q,1);
        idx = find(isnan(temp_mask)) ; %land mask: temp==0
        land_idx(idx) = 1 ;
        land_mask = find(land_idx ==1);
        ocean_mask = find(land_idx ==0);
        ocean_mask_name=strcat(workdir, '/test/input/', model_name, '_', 'ocean_mask_layer_', num2str(j,'%02i'), '.txt');
        save(ocean_mask_name,'ocean_mask','-ascii');

        %% extract temp data(exclude NaN)
        data =zeros(n*q- length(land_mask),m);
        for i = 1:m 
            temp2=squeeze(temp(:,:,i));
        %temp2(i,j,:,:)=permute(temp(:,:,j,i), [4 3 2 1]); %temp : lon, lat , depth, time
            temp2(land_mask) = [];
            data(:,i) = temp2(:)' ; % [space time]
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
        tt = strcat(workdir, '/test/input/', model_name, '_', var_name, '_', num2str(j,'%02i'), '.data');
        save(tt,'data2','-ASCII');
        disp(['dimension of sampling stations : ',num2str(slen)]);
        disp(['number of sampling points at each station : ',num2str(tlen)]);
    end



    %% STEP 2_1 Makescript_eigen

    %% make cseof directory for each variable, model
    cd([workdir, '/test/']);
    cseof_output_dir = [workdir, '/test/output'];
    system(['rm -rf ', cseof_output_dir, '/*']);
    mkdir(cseof_output_dir);

    name= cell(28);
    %----------STEP2. make script file-----------------------------------%
    name{1} = '#!/bin/csh';     %use csh shell
    name{2} = char('');
    name{3} = char(['gfortran  -fno-backtrace -o eigen ', workdir, '/programs/eigen/eigenx.f']); %compile the eigen file
    name{4} = char('');
    name{5} = char('cat >! eof.com <<ENDc'  );
%     ttt = ['./../data/', model_name, '_', var_name, '_', num2str(j,'%02i'), '.data'];
%     ttt = ['./', model_name, '_', var_name, '_', num2str(j,'%02i'), '.data'];
%     system(['ln -sf ', tt, ' ', ttt]);  %% eigen fortran codes cannot recognize filenames longer than 50 characters.
%     name{6} = char(ttt);
    name{6} = char(tt);  % input file name
    name{7} = char(['((', num2str(tlen), 'e16.7))'] );  %ascii file read, (time_num)e16.7
    name{8} = char([num2str(slen), ' 1']); %space number
    name{9} = char(num2str(tlen));  %time number
    name{10} = char('1            '    );  %space*time matrix or time*space matrix in matlab(1 or 2)
    name{11} = char('0            ' );                            
    name{12} = char('1          '  );
    name{13} = char('0              ' ); %% area adjustment (0: No)
    name{14} = char('15. 37.' );  %% starting latitude and increment for area adjustment
%     name{15} = char('99.99');  %want to modes to explain 95% (percent variance)
    name{15} = char('100.0');  %want to modes to explain 95% (percent variance)
    name{16} = char('1.' );  %% EOF scaling factor
    eofmodenum= 50;
    maxmodenum = num2str(eofmodenum);
    name{17} = char(maxmodenum);  %% number of EOFs to be printed                              
    name{18} = char('0');  %% PC normalization
    eof_tt = strcat(cseof_output_dir, '/eof_', model_name, '_', var_name, '.dat');
    name{19} = char(eof_tt);  %save eof LV file    
    name{20} = char('DIR');                                       %DIR : binary
    pct_tt = strcat(cseof_output_dir, '/pct_', model_name, '_', var_name, '.dat');
    name{21} = char(pct_tt);  %save eof pct file
    name{22} = char('(5e13.5)' );     %% ASCII format "in parenthesis"                               %save option: 6e13.5 acsii
    name{23} = strtrim(char('ENDc'));
    name{24} = char(['./eigen < ./eof.com'] );
    name{25} = char(['mv -f ', workdir, '/test/scripts/inform.d ', cseof_output_dir, '/inf_', model_name, '_', var_name,'.d']); %save information file
    name{26} = char(['mv -f ', workdir, '/test/scripts/avg.d ', cseof_output_dir, '/avg_', model_name, '_', var_name,'.d']);  %save average file
    name{27} = char('rm -f eigen eof.com' );
%     name{27} = char('rm -f ei' );
    name{28} = char('' );
    %     name = char(name,strtrim(name1));
    % end

    fid = fopen([workdir, '/test/scripts/eigen.c'], 'w+')
    for nline=1:length(name)
        fprintf(fid, '%s\n', name{nline});
    end
    fclose(fid);
%     system(['cd ', workdir, '/test/scripts'])
    cd([workdir, '/test/scripts']);
    system(['export PATH=$PATH:/usr/local/bin; csh -xv eigen.c > ', cseof_output_dir, '/eigen_log.log']);
        
    %% STEP 2_1 Makescript_cseof
    name = cell(30);
    name{1} = '#!/bin/csh';     %use csh shell
    % for i = 1:length(nn);
    name{2} = char('');
    name{3} = char( ['gfortran -fno-backtrace -o cseof ', workdir, '/programs/eigen/cseof.f']); %compile the eigen file
    name{4} = char( '' );
    name{5} = char( 'cat >! cseof.com <<ENDc'  );
    name{6} = char('0');  %% job number (only 0)
    name{7} = char( pct_tt  );     %load variable file
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
%     name{17} = char( '99.99' );   %want explain %
    name{17} = char( '100.0' );   %want explain %
    cseofmodenum=eofmodenum;
    name{18} = char( num2str(cseofmodenum) );    %want max.mode number
    name{19} = char( '1.');   %% eof scale factor (Do not change)
    name{20} = char( '2');  %% 1: rc ts    2 : cov
    name{21} = char( 'ENDc' );
    name{22} = char( [workdir, '/test/scripts/cseof < ', workdir, '/test/scripts/cseof.com'] );
    cinf_tt = strcat(cseof_output_dir, '/cinf_', model_name, '_', var_name, '.d');
    blo_tt = strcat(cseof_output_dir, '/blo_', model_name, '_', var_name, '.d');
    cpct_tt = strcat(cseof_output_dir, '/cpct_', model_name, '_', var_name, '.d');
    ceig_tt = strcat(cseof_output_dir, '/ceig_', model_name, '_', var_name, '.d');
    covm_tt = strcat(cseof_output_dir, '/covm_', model_name, '_', var_name, '.d');
    hcoef_tt = strcat(cseof_output_dir, '/hcoef_', model_name, '_', var_name, '.d');
    name{23} = char( ['mv -f inform.d ', cinf_tt]);
    name{24} = char(  ['mv -f Bloch.d ', blo_tt] );
    name{25} = char(  ['mv -f pcts.d ', cpct_tt] );
    name{26} = char(  ['mv -f eigen.d ', ceig_tt] );
    name{27} = char(  ['mv -f covm.d ', covm_tt] );
    name{28} = char(  ['mv -f hcoef.d ', hcoef_tt] );
    name{29} = char( 'rm -f cseof cseof.com' );
    name{30} = char( '' );
    % end

    fid = fopen([workdir, '/test/scripts/cseof.c'], 'w+')
    for nline=1:size(name,1)
        fprintf(fid, '%s\n', name{nline});
    end
    fclose(fid);
    system(['export PATH=$PATH:/usr/local/bin; csh -xv cseof.c > ', cseof_output_dir, '/cseof_log.log']);
    
    
    
    
    %% STEP 2_5 Makescript_recast
    %------case1. one layer-----------------%
    name = cell(23);
    name{1} = '#!/bin/csh';     %use csh shell

    name{2} = char('');
    name{3} = char(['gfortran -fno-backtrace -o recast ', workdir, '/programs/util/recastx.f']);
    name{4} = char( '' );
    name{5} = char( 'cat >! recast.com <<ENDc'  );
    name{6} = char( eof_tt  ); %predictor eof file
    name{7} = char( 'DIR' );
    name{8} = char( 'nofile'  );
    name{9} = char( num2str(cseofmodenum));                     %max. eof mode number
    name{10} = char(  num2str(slen) );   %space number
    name{11} = char(  '1.      ' );                              
    name{12} = char(  [num2str(slen), ' 1']  );     %space structure
    % name{13} = char( strcat('../regress/blo_reg_', 'NorESM1-M_', varname,'.d') ); %regressed blo file
    name{13} = char( blo_tt ); %blo file
    name{14} = char( 'SEQ' );
    name{15} = char( '0');
    sizeLV = LVnumber * cseofmodenum;
    name{16} = char( num2str(sizeLV) );  %LV number * cseof mode number (365*cseofmodenum)
    name{17} = char( '1' );
    LV_tt = strcat(cseof_output_dir, '/LV_', model_name, '_', var_name, '.dat');
    name{18} = char( LV_tt);     %LV file name
    name{19} = char( 'DIR');
    name{20} = char( 'ENDc' );
    name{21} = char( [workdir, '/test/scripts/recast < ', workdir, '/test/scripts/recast.com'] );
    name{22} = char( 'rm -f recast recast.com' );
    name{23} = char( '' );

    fid = fopen([workdir, '/test/scripts/recast.c'], 'w+')
    for nline=1:size(name,1)
        fprintf(fid, '%s\n', name{nline});
    end
    fclose(fid);
    system(['export PATH=$PATH:/usr/local/bin; csh -xv recast.c > ', cseof_output_dir, '/recast_log.log']);
    
    cd(presentdir);
    
    
    %% STEP 3 LV_layer_merge
%     filepath = [workdir, '/test/output'];
%     filename = [var_name, '_interp_', model_name, '_', scen_name, ...
%         '_', num2str(tempyear,'%04i'), '.nc']

%     file = strcat(filepath, '/', filename)

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
        filepath = [workdir, '/test/output']
        LV_name = strcat(filepath,'/', 'LV_', model_name, '_', var_name, '.dat');
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
    figure_output_dir = [workdir, '/test/figure'];
    mkdir([figure_output_dir]);
% % %     %load loading vector file
% % %     load([cseof_output_dir, '/lv_layer_merge.mat']);

    mode = eofmodenum;
    %load cseof pc time series
    pct_data = importdata(cpct_tt)';
    pct = pct_data(:)';
    pct = reshape(pct,tlen,length(pct)/tlen);
    %pct(:,2) = -pct(:,2);
    %plot pc time series

    for nyear = 1:length(inputyear)
        tempyear = inputyear(nyear);
        for month=1:12
            xData((12*(nyear-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]); %% time
        end
    end
    
    %% unnecessary every pct plot
%     clear i
%     for i = 1:mode
%         plot(xData, pct(:,i));
%         tt = [var_name, '-',num2str(i, '%02i'),'mode'];
%         title(tt,'fontsize',18,'fontweight','bold');  
%         datetick('x', 'yy', 'keepticks')
%         xlabel('Time(year)','FontSize',20) ;
%         set(gca,'FontSize',14);
%         saveas(gcf, strcat([figure_output_dir, '/', tt, '.png']),'png');
%         close all
%         i
%     end

    LV_increment1(LV_increment1==0) = NaN;




    %%
    %% blue-white-red colormap
    %%
    
    i=1:20;
      bwrmap(i,1)= 0.;
      bwrmap(i,2)= 0.2:(0.3/19.):0.5;
      bwrmap(i,3)= 0.4:(0.6/19.):1.;
    
    i=21:49;
      bwrmap(i,1)= 0.:(1./28.):1.;
      bwrmap(i,2)= 0.5:(0.5/28.):1.;
      bwrmap(i,3)= 1.;   
    
    i=49:51;
      bwrmap(i,1)= 1.;
      bwrmap(i,2)= 1.;
      bwrmap(i,3)= 1.;
    
    i=51:56;
      bwrmap(i,1)= 1.;
      bwrmap(i,2)= 1.:(-0.1/5.):0.9;
      bwrmap(i,3)= 1.:(-0.1/5.):0.9;
    
    i=56:70;
      bwrmap(i,1)= 1.;
      bwrmap(i,2)= 0.9:(-0.45/14.):0.45;
      bwrmap(i,3)= 0.9:(-0.45/14.):0.45;
     
    i=70:80;
      bwrmap(i,1)= 1.;
      bwrmap(i,2)= 0.45:(-0.45/10.):0.;
      bwrmap(i,3)= 0.45:(-0.45/10.):0.;
     
    i=80:100;
      bwrmap(i,1)= 1.:(-0.6/20.):0.4;
      bwrmap(i,2)= 0.;
      bwrmap(i,3)= 0.;
      

    save('/Users/kimyy/CSEOF/data.mat');
    
    
    
    load('/Users/kimyy/CSEOF/data.mat');
    %% data reconstruction
    temp_recon=zeros(size(temp));
    for mi=1:12
        for modei=1:eofmodenum
            tmp_lv=repmat(squeeze(LV_increment1(:,:,:,mi,modei)), [1,1,length(inputyear)]);
            tmp_pct=reshape(pct(mi:12:size(pct,1)-12+mi,modei), [1,1,length(inputyear)]);
            temp_recon(:,:,mi:12:size(pct,1)-12+mi)=temp_recon(:,:,mi:12:size(pct,1)-12+mi) + ...
                tmp_lv .* tmp_pct;
        end
        temp_recon1(:,:,mi:12:size(pct,1)-12+mi)=temp_recon(:,:,mi:12:size(pct,1)-12+mi) + ...
            mean(temp(:,:,mi:12:size(pct,1)-12+mi),3)-mean(temp_recon(:,:,mi:12:size(pct,1)-12+mi),3);
    %     temp_recon2(:,:,mi:12:size(pct,1)-12+mi)=temp_recon(:,:,mi:12:size(pct,1)-12+mi) + ...
    %         mean(temp(:,:,:)-temp_recon(:,:,:),3);
    end

    % pcolor(temp(:,:,14)'-temp_recon(:,:,14)'); shading flat; colorbar;
    
    %practice
    %raw plot
    close all;
    hold off
    plot(squeeze(temp_recon1(45,35,1:60)))
    hold on
    plot(squeeze(temp(45,35,1:60)))
    
    % get climatological mean
    temp_clim=reshape(temp,[45,39,12,40]);
    temp_clim_cycle=repmat(mean(temp_clim,4), [1,1,40]);
    % anomaly plot
    close all;
    hold off
    plot(squeeze(temp_recon1(45,35,:)-temp_clim_cycle(45,35,:)), 'linewidth', 2)
    hold on
    plot(squeeze(temp(45,35,:)-temp_clim_cycle(45,35,:)), 'linewidth', 2)
    
    % quantitative error check
    rmse_v=squeeze(sqrt(sum((temp_recon1-temp).^2,3)./480)); 
    std_ano=std(temp-temp_clim_cycle,1,3);
    std_raw=std(temp,1,3);
    std_recon_ano=std(temp_recon1-temp_clim_cycle,1,3);
    
    
    %% main fig 1
    close all;
    sb1=subplot(4,2,1); % mean of raw
    pcolor(mean(temp,3)'); shading flat; colorbar; colormap(sb1,jet); caxis([5 28]); set(gca, 'fontsize', 20);
    sb2=subplot(4,2,2); % mean of reconstructed data
    pcolor(mean(temp_recon1,3)'); shading flat; colorbar; colormap(sb2,jet); caxis([5 28]); set(gca, 'fontsize', 20);
    sb3=subplot(4,2,3); % mean bias
    pcolor(mean(temp-temp_recon1,3)'); shading flat; colorbar; colormap(sb3,bwrmap); set(gca, 'fontsize', 20);
    sb4=subplot(4,2,4); % error at certain time
    pcolor(temp(:,:,115)'-temp_recon1(:,:,115)'); shading flat; colorbar; colormap(sb4,bwrmap); set(gca, 'fontsize', 20);
    sb5=subplot(4,2,5); % raw time series at certain point
    plot(squeeze(temp_recon1(45,35,1:60))); hold on; plot(squeeze(temp(45,35,1:60)));
    sb6=subplot(4,2,6); % anomaly time series at certain point
    plot(squeeze(temp_recon1(45,35,:)-temp_clim_cycle(45,35,:)), 'linewidth', 2); hold on; plot(squeeze(temp(45,35,:)-temp_clim_cycle(45,35,:)), 'linewidth', 2)
    sb7=subplot(4,2,7); % rate of std between recon and raw
    pcolor((std_recon_ano./std_ano)'); shading flat; colorbar; colormap(sb7,flip(autumn)); set(gca, 'fontsize', 20);
    sb8=subplot(4,2,8); % rmse
    pcolor(rmse_v'); shading flat; colorbar; colormap(sb8,flip(spring)); set(gca, 'fontsize', 20);
    
    % pcolor(sqrt(temp-temp_recon1,3)'); shading flat; colorbar; colormap(sb3,bwrmap); set(gca, 'fontsize', 20);
    
    
    
    
    
    %% main fig 2
    close all;
    sb1=subplot(2,1,1);
    plot(pct(:,1));
    sb2=subplot(2,1,2);
    plot(pct(:,2));



end  % for var_name