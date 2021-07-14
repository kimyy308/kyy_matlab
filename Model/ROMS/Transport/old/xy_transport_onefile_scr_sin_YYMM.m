% function [total_volume,total_salt,total_fresh] = xy_transport_onefile(filename,g,selectdepth);
% ====================================================================
% [total_volume,total_salt,total_fresh]= xy_transport_onefile(filename,g,selectdepth);
% calculate transports of  volume(in Sv, 10^6 m^3/sec)
%                          heat  (in PW, 10^15 W)
%                          salt  (in 10^9 kg)
%
% across a line segment (slice along a constant I or J)
%
% 'filename'       = history or average file name
% 'grid_nick_name' = grid information such as 'eas' 'hudson' 'latte'
% 'selectdepth' = surface to which depth (m) such as 100, 500, 50000
%
% keep in mind that
% vertical coordinate changes in time = h + zeta(t) in ROMS
%
% USE: xy_transport_function.m
%
% ====================================================================
% onefile   for single   input file  and single segment.
% BJ Choi, Marhch07, 2006.
clear all; clc; tic;
%   if( nargin < 1)
%     filename='./../avg_his/avg_eRTD_0002.nc';
%     filename='./../monthly/avg_Ryw10nao_notc_1993_01.nc'; 
    filename='C:\Users\kyy\Desktop\nwp\TRANSP~1\ocean_avg_0409.nc';     

% year_start = 1993 ;
% year_end   = 2005 ;
% % year_end   = 1998 ;
% mon_number = [1 12]; %% the last year month number

year_start = 1992 ;
year_end   = 2008 ;
year_start = 1980 ;
year_end   = 2009 ;

year_start = 2003 ;
year_end   = 2011 ;
mon_number = [1 12 ]; %% the last year month number

header='';
% header='ref35_';
% header='ref34p7_';
% header='D2000m_';

%     disp(['Open file name = ',filename])
        sfile=filename;
%     g = grd('yw_n10_4m');
%     g = grd('yw_skku_curve_0530');
%     g = grd('yw_skku_curve_1001');
    g=grd('NWP');
    
    selectdepth= 200;
    selectdepth= 2000;    
%     selectdepth  = input('surface to which depth (m) = ');
%   end

    s_max = 32.76; % maximum salinity for the calculation of 
                   % freshwater flux - Hudson River, NY Bight.

    s_max = 34.00; % maximum salinity for the calculation of 
                   % freshwater flux - Hudson River, NY Bight.
    s_max = 34.96522443559434;
    s_max = 34.8666786780724;
    s_max = 34.57;
    s_max = 34.5629306519175 ; % spinup
    s_max = 34.5586083395882 ; % real ~ 103
%     s_max = 34.5535161339153
    s_max = 34.5581101007304 ;%yw10 - 9302
    s_max = 34.5771900087778 ;%yw10_re - 9305
    s_max = 34.5766078608146 ; % yw10 - notc
    s_max = 34.5765824422818 ; %notc_nao
    s_max = 34.571722026517 ; % notc_nao_skku1 
 s_max = 34.5730282399533   ; % notc_nao_skku_re2
  s_max = 34.573    ; % notc_nao_skku_re2
%   s_max = 34.65   ; % notc_nao_skku_re2  
%   s_max = 34.57   ; % notc_nao_skku_re2    

stnum = 1;
stnum = (year_start-1980)*12 +1;
% *************************************************************
%
%   END OF USER DEFINED VARIABLES
%
% *************************************************************

  % size of grids
  [r,c] = size ( g.lon_rho );
  mask3d_rho=repmat(g.mask_rho,[1 1 g.N]);
  mask3d_rho=permute(mask3d_rho,[3 1 2]);

  [r1,c1] = size ( g.lon_u );
  mask3d_u=repmat(g.mask_u,[1 1 g.N]);
  mask3d_u=permute(mask3d_u,[3 1 2]);

  [r2,c2] = size ( g.lon_v );
  mask3d_v=repmat(g.mask_v,[1 1 g.N]);
  mask3d_v=permute(mask3d_v,[3 1 2]);
  
  % transport from surface to which depth (m)

    if ( selectdepth > 0 ) 
     selectdepth = selectdepth*-1;
    end

t_point=3;  % kor, tsu, soy, tiw_e, tiw_w
point_name={'Korea Strait','YS','CJ','Taiwan Strait','Kuro intrusion'};
%     'Tsugaru Strait','West Taiwan Strait','East Taiwan Strait'};
% kts = [128.5696  35.185 130.1933  33.2432   ...
% kts = [128.8803  35.5872 129.9648  33.1279   ...
% kts = [     128.6 34.85  129.5  33.2265  ... KS
kts = [  128.79  34.8999  129.7  33.3353  ... KS..e.w
    119.6354  34.5195 126.4802  34.5195   ...   YS
    120.3889  32.3077  120.1667  31.742  ...   CJ   
    119.3894  25.9829  121.4576  25.1472  ...
    129.5  33.2265  121.658  24.9367  ];

%     120.3889  32.3077  120.1667  31.7422  ...   CJ  
%     120.9444  32.2137  120.8333  31.5529  ...   CJ  -- OK
        
% kts  = [ 128.25 35.25 132.75 35.25 ...  % KS     
%          117.25 24.25 121.25 24.25 ...  % TS     
%          121.25 24.25 126.25 24.25 ...  % Kuro   
%          ];


%        118.5954  25   120.9149  25    121.3789 25    125.0902 25]; %mks 2

    %choose two grid points (lon, lat)
    
%     endpt_lon = [128.5696 130.1933 142.2552 140.8634 ];
%     endpt_lat = [ 35.185   33.2432  46.4677  40.2885 ];

% [M N] =size(endpt_lon);
% for point=1:N/2
 t_point = length(kts)/4 ;
 
figure
pcolorjw(g.lon_rho,g.lat_rho,g.mask_rho_nan)
for k =1: t_point
    kk = 4*(k-1) ;
    line([kts(kk+1) kts(kk+3)],[kts(kk+2) kts(kk+4)],'Marker','.','LineStyle','-')
end
saveas(gcf,'connect_3point','png')   

% aaaa
 
depth = selectdepth;
sn=0;

for i = 1:20   %%%%
masku(i,:,:) = g.mask_u ;
maskv(i,:,:) = g.mask_v ;
end

for st=1:t_point
     
    disp([char(point_name(st)),' Depth to ',num2str(depth)])
     endpt_lon(1) = kts(sn+1); 
     endpt_lat(1) = kts(sn+2);
     endpt_lon(2) = kts(sn+3);
     endpt_lat(2) = kts(sn+4);
     sn=sn+4;
    
  for cpts=1:2 % corner points
    dist = sqrt(  ( g.lon_rho - endpt_lon(cpts) ).*( g.lon_rho - endpt_lon(cpts) ) + ...
                  ( g.lat_rho - endpt_lat(cpts) ).*( g.lat_rho - endpt_lat(cpts) ) );
    ind=find( min( dist(:) ) == dist );
    % closest points row and column indice
    row_index = mod ( ind - 1, r ) + 1;
    col_index = floor( (ind - 1) / r ) + 1;
    corner_endpt_col(cpts)=col_index;
    corner_endpt_row(cpts)=row_index;
  end


  % my xy_transport_onefile works only if corner_endpt_row(2) < corner_endpt_row(1).
   if( corner_endpt_row(2) > corner_endpt_row(1)  )
     tmp_col=corner_endpt_col(2);
     tmp_row=corner_endpt_row(2);
     corner_endpt_col(2)=corner_endpt_col(1);
     corner_endpt_row(2)=corner_endpt_row(1);
     corner_endpt_col(1)=tmp_col;
     corner_endpt_row(1)=tmp_row;
     beep
     disp(' === switching two end points === ')
   end

  % longitude and latitude coordinate.

    for i=1:length(endpt_lat)
     xx(i)=g.lon_rho(corner_endpt_row(i),corner_endpt_col(i));
     yy(i)=g.lat_rho(corner_endpt_row(i),corner_endpt_col(i));
    end
    distance_r = m_lldist ( xx, yy );

%     figure(f)
%     hold on
%     plot(xx,yy,'x-k','LineWidth',2)
%     colorbar
%     title('bottom topography (m)')

%  transect information
 
    % delj = j increasment
    if( corner_endpt_col(2) >= corner_endpt_col(1) ) 
    delj=1;
    west2east_transect=1; % previously zonaltransect
    else
    delj=-1;
    west2east_transect=0; % previously meridionaltransect
    end

    % deli = i increasment
    if( corner_endpt_row(2) > corner_endpt_row(1) ) 
    deli=1;
    else
    deli=-1;
    end


    %               i j
    %             row col
    % g.lon_rho: [142x254 double] for latte grid
    % g.lon_u:   [142x253 double]
    % g.lon_v:   [141x254 double]

    xzero=g.lon_rho( corner_endpt_row(1), corner_endpt_col(1) );
    yzero=g.lat_rho( corner_endpt_row(1), corner_endpt_col(1) );
    xone=g.lon_rho( corner_endpt_row(2), corner_endpt_col(2) );
    yone=g.lat_rho( corner_endpt_row(2), corner_endpt_col(2) );
    slope=( yone-yzero) / (xone - xzero);
    % A x + B y + C = 0;
    A=slope;
    B=-1;
    C=-slope*xzero+yzero;
    D=sqrt( A*A + B*B );
    % distance = abs( A x + B y + C ) / D

%   grid information

   %N  is the number of vertical levels
   %hz is thickness  of each level
    N = g.N;  
    [M L]=size(g.h);
    hz=g.z_w(2:N+1,:,:)-g.z_w(1:N,:,:); % z_w: [31x142x254]
    dx = 1./g.pm;
    dy = 1./g.pn;
    dx_v=0.5*(dx(1:M-1,:)+dx(2:M,:));
    dy_u=0.5*(dy(:,1:L-1)+dy(:,2:L));
    g_h_v=0.5*(g.h(1:M-1,:)+g.h(2:M,:));
    g_h_u=0.5*(g.h(:,1:L-1)+g.h(:,2:L));

    avgdxdy = mean([ mean( mean( dx ) ), mean( mean( dy ) ) ]);

fnum=0;
%   load model output file
%      for fileloop=file_start:file_step:file_end
%         ll=length(num2str(fileloop));
%         filename(end-2-ll:end-3)=num2str(fileloop); 
        
for yy = year_start : year_end ;
    filename(end-9:end-6)=num2str(yy);
mon_start = 1  ; %% the default year month number
mon_end   = 12 ; %% the default year month number 

    if yy == year_end
        mon_start = mon_number(1);
        mon_end   = mon_number(end);
    end
for mm =  mon_start : mon_end
    filename(end-4:end-3)=num2char(mm,2); 
        fnum=fnum+1;
    nc=netcdf(filename,'read');
%     disp([' opening your data file: ', filename])
    zeta=nc{'zeta'}(:); % zeta(time, eta_rho, xi_rho)
    u=nc{'u'}(:);       % u(time, s_rho, eta_u, xi_u)
    v=nc{'v'}(:);       % v(time, s_rho, eta_u, xi_u)
                u = u.*masku ; v = v.*maskv ;   %%%%
    temp=nc{'temp'}(:); % temp(time, s_rho, eta_rho, xi_rho)
    salt=nc{'salt'}(:); % salt(time, s_rho, eta_rho, xi_rho)
    close(nc)

    salt = salt.*mask3d_rho; % zero for land, ** very important **

    szero = ones(size(salt)).*s_max;
    fresh = ( szero - salt ) ./ szero; % freshwater fraction
    fresh = fresh.*mask3d_rho; % zero for land, ** very important **

%   vertical coordinate changes in time 
%   because sea surface height changes in time.
%   thickness of each layer changes propotional to total water thickness.

    h_total = g.h + zeta;       %total water thickness
    for level=1:N               %thickness of each layer
     Hz(level,:,:)=squeeze(hz(level,:,:)).*(h_total./g.h);
    end

    % average Hz to  Arakawa-C u points 

    Hz_u=0.5*(Hz(:,:,1:L-1)+Hz(:,:,2:L)); % each layer thickness
     z_u(1,:,:)=-g_h_u(:,:);             % z @ bottom of each layer
    for k=2:+1:N
     z_u(k,:,:)=z_u(k-1,:,:)+Hz_u(k-1,:,:);
    end

    temp_u=0.5*(temp(:,:,1:L-1)+temp(:,:,2:L)); % each layer temp at u point
    salt_u=0.5*(salt(:,:,1:L-1)+salt(:,:,2:L)); % each layer salt at u point
    fresh_u=0.5*(fresh(:,:,1:L-1)+fresh(:,:,2:L)); % each layer freshwater at u point

    % average Hz to  Arakawa-C v points 

    Hz_v=0.5*(Hz(:,1:M-1,:)+Hz(:,2:M,:)); % each layer thickness
     z_v(1,:,:)=-g_h_v(:,:);             % z @ bottom of each layer
    for k=2:+1:N
     z_v(k,:,:)=z_v(k-1,:,:)+Hz_v(k-1,:,:);
    end

    temp_v=0.5*(temp(:,1:M-1,:)+temp(:,2:M,:)); % each layer temp at u point
    salt_v=0.5*(salt(:,1:M-1,:)+salt(:,2:M,:)); % each layer salt at u point
    fresh_v=0.5*(fresh(:,1:M-1,:)+fresh(:,2:M,:)); % each layer freshwater at u point


%   ====================================================================================
%   find path from corner_endpt(1) to corner_endpt(2)

    icount=1;
       col_index(icount)=corner_endpt_col(1); 
       row_index(icount)=corner_endpt_row(1);
       on_vpoint(icount)=0;
       vpoint=0;
       xpoint=g.lon_u( row_index(icount), col_index(icount) );
       ypoint=g.lat_u( row_index(icount), col_index(icount) );
     
    signline(icount)=(ypoint - yzero) - slope*(xpoint -xzero);
    dist(icount)= abs( A*xpoint + B*ypoint + C ) / D;
    tmp_dist=dist(icount);
    dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone]);
    flag_approach = 1;

    %  icount
    %  row_index(icount)
    %  col_index(icount)
    %  vpoint
    %  signline(icount)
    %  dist2endpoint(icount)
%       plot(  xpoint,  ypoint , 'ro')


  if (west2east_transect)

    while ( dist2endpoint(icount) > avgdxdy  &&  flag_approach )

     icount=icount+1;

     if ( vpoint == 1 )

       col_index(icount)=col_index(icount-1)+delj;
       if ( on_vpoint(icount-1) == 1)
        row_index(icount)=row_index(icount-1);
       else
        row_index(icount)=row_index(icount-1)+deli;
       end
       xpoint=g.lon_v( row_index(icount), col_index(icount) );
       ypoint=g.lat_v( row_index(icount), col_index(icount) );
       signline(icount)=(ypoint - yzero) - slope*(xpoint -xzero);
       dist(icount)= abs( A*xpoint + B*ypoint + C ) / D;
    
       if ( signline(icount)*signline(icount-1) < 0  ...
            ||   dist(icount) <= dist(icount-1)       ...
            ||   dist(icount) <= tmp_dist                  )
         tmp_dist=0;
         on_vpoint(icount)=1;
         plot(  xpoint,  ypoint , 'ro')
         dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone]);
       else
         tmp_dist=dist(icount);
         vpoint=0;
         icount=icount-1;
         plot(  xpoint,  ypoint , 'gx')
       end

      else % on upoint

       col_index(icount)=col_index(icount-1);
       if ( on_vpoint(icount-1) == 0)
       row_index(icount)=row_index(icount-1)+deli;
       else
       row_index(icount)=row_index(icount-1);
       end
       xpoint=g.lon_u( row_index(icount), col_index(icount) );
       ypoint=g.lat_u( row_index(icount), col_index(icount) );
       signline(icount)=(ypoint - yzero) - slope*(xpoint -xzero);
       dist(icount)= abs( A*xpoint + B*ypoint + C ) / D;
    
       if (      signline(icount)*signline(icount-1) < 0 ...
            ||   dist(icount) <= dist(icount-1)   ...
            ||   dist(icount) <= tmp_dist                       )
         tmp_dist=0;
         on_vpoint(icount)=0;
         plot(  xpoint,  ypoint , 'ro')
         dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone]);
       else
         tmp_dist=dist(icount);
         vpoint=1;
         icount=icount-1;
         plot(  xpoint,  ypoint , 'gx')
       end

      end % if ( on_vpoint == 1 )

      if( icount > 3 &&  dist2endpoint(icount) > dist2endpoint(icount-3) )
       flag_approach = 0;
      end

     % icount
     % row_index(icount)
     % col_index(icount)
     % vpoint
     % signline(icount)
     % dist2endpoint(icount)

    end % while

  else % if (west2east_transect)

    while ( dist2endpoint(icount) > avgdxdy  &&  flag_approach )

     icount=icount+1;

     if ( vpoint == 1 )

       if ( on_vpoint(icount-1) == 1)
        col_index(icount)=col_index(icount-1)+delj;
        row_index(icount)=row_index(icount-1);
       else
        col_index(icount)=col_index(icount-1);
        row_index(icount)=row_index(icount-1)+deli;
       end

       xpoint=g.lon_v( row_index(icount), col_index(icount) );
       ypoint=g.lat_v( row_index(icount), col_index(icount) );
       signline(icount)=(ypoint - yzero) - slope*(xpoint -xzero);
       dist(icount)= abs( A*xpoint + B*ypoint + C ) / D;
    
       if ( signline(icount)*signline(icount-1) < 0  ...
            ||   dist(icount) <= dist(icount-1)       ...
            ||   dist(icount) <= tmp_dist                  )
         tmp_dist=0;
         on_vpoint(icount)=1;
         plot(  xpoint,  ypoint , 'ro')
         dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone]);
       else
         tmp_dist=dist(icount);
         vpoint=0;
         icount=icount-1;
         plot(  xpoint,  ypoint , 'gx')
       end

      else % on upoint

       if ( on_vpoint(icount-1) == 0)
        col_index(icount)=col_index(icount-1);
       row_index(icount)=row_index(icount-1)+deli;
       else
        col_index(icount)=col_index(icount-1)+delj;
       row_index(icount)=row_index(icount-1);
       end


       xpoint=g.lon_u( row_index(icount), col_index(icount) );
       ypoint=g.lat_u( row_index(icount), col_index(icount) );
       signline(icount)=(ypoint - yzero) - slope*(xpoint -xzero);
       dist(icount)= abs( A*xpoint + B*ypoint + C ) / D;
    
       if (      signline(icount)*signline(icount-1) < 0 ...
            ||   dist(icount) <= dist(icount-1)   ...
            ||   dist(icount) <= tmp_dist                       )
         tmp_dist=0;
         on_vpoint(icount)=0;
         plot(  xpoint,  ypoint , 'ro')
         dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone]);
       else
         tmp_dist=dist(icount);
         vpoint=1;
         icount=icount-1;
         plot(  xpoint,  ypoint , 'gx')
       end

      end % if ( on_vpoint == 1 )

      if( icount > 3 &&  dist2endpoint(icount) > dist2endpoint(icount-3) )
       flag_approach = 0;
      end

     %icount
     % row_index(icount)
     % col_index(icount)
     % vpoint
     % signline(icount)
     % dist2endpoint(icount)

    end % while

   end % if (west2east_transect)


    total_volume=0;
    total_salt=0;
    total_fresh=0;
 
    for index=1:icount

      vpoint=on_vpoint(index);    
      xy_transport_function     
      
      if( west2east_transect == 0 && vpoint == 1 )
        total_volume = total_volume - sum_segment;
        total_salt   = total_salt   - sum_segment_salt;
        total_fresh  = total_fresh  - sum_segment_fresh;
        freshtransect(index)=-sum_segment_fresh;
        voltransect(index)=-sum_segment;
      else
        total_volume = total_volume + sum_segment;
        total_salt   = total_salt   + sum_segment_salt;
        total_fresh  = total_fresh  + sum_segment_fresh;
        freshtransect(index)=sum_segment_fresh;
        voltransect(index)=sum_segment;
      end


    end
    disp(['YY/MM = ',num2str(yy),'/',num2char(mm,2),' transport ',char(point_name(st)),': ',num2str(total_volume/1e+6)])  
    trans(fnum,st) = total_volume/1e+6;
    salt_tr(fnum,st)  = total_salt/1e+6;
    fresh_tr(fnum,st)  = (total_fresh/1e+6)*1000;
     end   % for mm
end %% for yy
      filename=sfile;
end  % for st=1:t_point


% save -ascii transport_new_1000m.dat trans
% save -ascii salt_trans_new_1000m.dat salt_tr
% save -ascii fresh_trans_new_1000m.dat fresh_tr

idx = [stnum:stnum+fnum-1]  ;

trans1    = [ idx'   trans    ];
salt_tr1  = [ idx'   salt_tr  ];
fresh_tr1 = [ idx'   fresh_tr ];

outfile = [header,'1transport_connect.txt'];
disp([outfile])
fid = fopen(outfile,'w+');
fprintf(fid,'%%time     KS     YS    CJ     West_Taiwan  Kuro_pre \n');
fprintf(fid,'%4d %10.4f %10.4f %10.4f %10.4f %10.4f \n',trans1');
fclose(fid);

outfile = [header,'1trans_salt_connect.txt'];
disp([outfile])
fid = fopen(outfile,'w+');
fprintf(fid,'%%time     KS     YS    CJ     West_Taiwan  Kuro_pre \n');
fprintf(fid,'%4d %10.4f %10.4f %10.4f %10.4f %10.4f \n',salt_tr1');
fclose(fid);

outfile = [header,'1trans_fresh_connect.txt'];
disp([outfile])
fid = fopen(outfile,'w+');
fprintf(fid,'%%time     KS     YS    CJ     West_Taiwan  Kuro_pre \n');
fprintf(fid,'%4d %10.4f %10.4f %10.4f %10.4f %10.4f \n',fresh_tr1');
fclose(fid);



cal_mean_salt = 1 ;
cal_mean_salt = 0 ;

if cal_mean_salt == 1 
outfile = ['../1transport_connect.txt'];
V = load(outfile); idx = V(:,1)'; trans = V(:,2:end);
outfile = ['../1trans_salt_connect.txt'];
V = load(outfile); idx = V(:,1)'; salt_tr = V(:,2:end);
outfile = ['../1trans_fresh_connect.txt'];
V = load(outfile); idx = V(:,1)'; fresh_tr = V(:,2:end);

salt_tr2  = [ idx'   salt_tr./trans  ];
fresh_tr2 = [ idx'   fresh_tr./trans ];

outfile = ['../1trans_mean_salt_connect.txt'];
disp([outfile])
fid = fopen(outfile,'w+');
fprintf(fid,'%%time     KS     YS    CJ     West_Taiwan  Kuro_pre \n');
fprintf(fid,'%4d %10.4f %10.4f %10.4f %10.4f %10.4f \n',salt_tr2');
fclose(fid);

outfile = ['../1trans_mean_fresh_connect.txt'];
disp([outfile])
fid = fopen(outfile,'w+');
fprintf(fid,'%%time     KS     YS    CJ     West_Taiwan  Kuro_pre \n');
fprintf(fid,'%4d %10.4f %10.4f %10.4f %10.4f %10.4f \n',fresh_tr2');
fclose(fid);
end





% DD = [ [1:fnum]' trans ];
% outfile = '1_transport.dat';
% disp(['  writing  = ',outfile])
% fid = fopen(outfile,'w+');
% fprintf(fid,'%ic    KS   TS    Kuro \n');
% fprintf(fid,'%5d %8.3f %8.3f %8.3f \n',DD');
% fclose(fid);
% 
% DD = [ [1:fnum]'+length(1980:1991)*12  trans ];
% outfile = '1_transport_1980.dat';
% disp(['  writing  = ',outfile])
% fid = fopen(outfile,'w+');
% fprintf(fid,'%ic    KS   TS    Kuro \n');
% fprintf(fid,'%5d %8.3f %8.3f %8.3f \n',DD');
% fclose(fid);





%   if( west2east_transect )
%    disp(' from West to East transect ')
%   else
%    disp(' from East to West (reverse) transect ')
%   end
%    disp(' ')
% 
%   disp([' volume transport = ', num2str(total_volume) ,' m^3/s '])
%   disp(['                  = ', num2str(total_volume/1e+6) ,' x 10^6m^3/s '])
% 
%   disp([' salt volume transport = ', num2str(total_salt) ,' g/s '])
%   disp(['                  = ', num2str(total_salt) ,' x Kg/s '])
% 
%   disp([' freshwater volume transport = ', num2str(total_fresh) ,' m^3/s '])
%   disp(['                  = ', num2str(total_fresh/1e+6) ,' x 10^6m^3/s '])
% 
%   disp([' transect distance = ', num2str( distance_r/1000 ),' km '])
%   disp([' first  point, lon, lat = ',  num2str( xx(1) ),'  ',num2str( yy(1) )])
%   disp([' second point, lon, lat = ',  num2str( xx(2) ),'  ',num2str( yy(2) )])
% 
% 
% 
%   iplot=1;
%   if (iplot)
%    figure
%    subplot(3,1,1)
%    plot([1:icount],voltransect/100,'g-x')
%    xlabel('index')
%    ylabel('volume trans')
%    grid
%    subplot(3,1,2)
%    plot([1:icount],freshtransect/100,'g-x')
%    xlabel('index')
%    ylabel('freshwater trans')
%    grid
%    subplot(3,1,3)
%    plot([1:icount],on_vpoint,'ro')
%    xlabel('index')
%    legend('1 for vpoint, 0 for upoint')
%    grid
%   end

toc