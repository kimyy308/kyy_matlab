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
% Updated 13-Oct-2020 by Yong-Yub Kim.
% ====================================================================

clc;close all;clear all;
warning off;

system_name=computer;
if (strcmp(system_name,'PCWIN64'))
    % % for windows
    dropboxpath='C:\Users\user\Dropbox'; %% SNU_desktop
    addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
    addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
%     addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Analysis\Figure\seo_nwp_1_10\run']));
elseif (strcmp(system_name,'GLNXA64'))
    % % for linux
%     dropboxpath='/home01/kimyy/Dropbox'; %% ROMS server
    dropboxpath='/home/kimyy/Dropbox'; %% DAMO server
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
%     addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Analysis/Figure/seo_nwp_1_10/run']));
end

experiment='nwp_1_20';
% all_testname2 = {'test2107', 'test2108', 'test2109'};
all_testname2 = {'test2117', 'test2118', 'test2119', 'test2120', 'test2121'};

% all_testname2 = {'test2109'};

inputyear = [1985:2014];

t_point= 8;  
point_name={'korea','tsugaru','soya','taiwan','kuro_intrusion','yellowsea', 'west', 'east'};
kts = [128.0  35.2   134.0  35.2   ...  
       140.5  42.0   140.54 41.0   ...
       142.0  47.0   142.04 45.0   ...
       118.65 24.6   120.2  23.75  ...
       121.95 25.1   130.15 31.35  ... 
       121.8  31.9   126.5  34.35  ...
       128.75 35.16  129.35 34.58  ...
       129.35 34.58  130.4  33.53];
lgd_point = ['%%korea  tsugaru   soya   taiwan  kuro_intru  yellowsea west   east\n'];

for testnameind2=1:length(all_testname2)

    testname=all_testname2{testnameind2};
%     switch testname
%         case {'test61', 'test62', 'test63', 'test64'}
%             drivename='H:\';
%         case {'test57', 'test58', 'test59', 'test60'}
%             drivename='I:\';
%         case {'test65', 'test66', 'test67', 'test68'}
%             drivename='G:\';
%         case {'ens08', 'ens09', 'ens10'}
%             drivename='E:\';
%     end
%     datadir = [drivename, 'Data\Model\ROMS\', experiment, '\'];
    
    datadir=['/data2/RCM/CMIP6/Model/ROMS/nwp_1_20/output/monthly/'];
    
    for yearind= 1:length(inputyear)
        year=inputyear(yearind);
%         inputdir = [datadir, testname, '\run\', num2str(year, '%04i'), '\']; %% for Desktop
        inputdir = [datadir, testname, filesep, 'run', filesep, 'packed_monthly', filesep, num2str(year, '%04i'), filesep]; %% for Desktop

        outputdir = [datadir, testname, filesep, 'run', filesep, 'transport_barot', filesep];
        if (exist(strcat(outputdir) , 'dir') ~= 7)
            mkdir(strcat(outputdir));
        end 
        
        outfile = strcat(outputdir, experiment, '_monthly_',testname,'_',num2str(year,'%04i'),'.txt');

        fid = fopen(outfile,'w+');
        fprintf(fid,lgd_point);
        fclose(fid);
        for month=1:12
            yearstr = num2str(year,'%04i');
            monthstr = num2str(month,'%02i');
            filename = strcat(inputdir, 'pck_', testname,'_monthly_',yearstr,'_',monthstr,'.nc');
    %         g = grd('NWP_1_20', testname);
            if (exist('g.lon_u' , 'var') ~= 1)  % [z lat lon]
                g.lon_rho=ncread(filename, 'lon_rho')';
                g.lat_rho=ncread(filename, 'lat_rho')';
                g.lon_u=ncread(filename, 'lon_u')';
                g.lat_u=ncread(filename, 'lat_u')';
                g.lon_v=ncread(filename, 'lon_v')';
                g.lat_v=ncread(filename, 'lat_v')';
                g.mask_rho=ncread(filename, 'mask_rho')';
                g.mask_u=ncread(filename, 'mask_u')';
                g.mask_v=ncread(filename, 'mask_v')';
                g.h=ncread(filename, 'h')';
                g.pm=ncread(filename, 'pm')';
                g.pn=ncread(filename, 'pn')';

                [len_lat_rho, len_lon_rho] = size ( g.lon_rho );
                g.dx = 1./g.pm;  % pm : curvilinear coordinate metric in XI
                g.dy = 1./g.pn;  % pn : curvilinear coordinate metric in ETA
                g.dx_v=0.5*(g.dx(1:len_lat_rho-1,:)+g.dx(2:len_lat_rho,:));
                g.dy_u=0.5*(g.dy(:,1:len_lon_rho-1)+g.dy(:,2:len_lon_rho));
                g.h_u=0.5*(g.h(1:len_lat_rho-1,:)+g.h(2:len_lat_rho,:));
                g.h_v=0.5*(g.h(:,1:len_lon_rho-1)+g.h(:,2:len_lon_rho));
    %             [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
    %             [lon_u_min, lon_u_max, lat_u_min, lat_u_max] = findind_Y(1/20, lonlat(1:4), lon_u, lat_u);
    %             [lon_v_min, lon_v_max, lat_v_min, lat_v_max] = findind_Y(1/20, lonlat(1:4), lon_v, lat_v);
    %             cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
    %             cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
            end
            zeta=ncread(filename, 'zeta')';
            zeta(zeta<=-1000)=NaN;
            zeta(zeta>=1000)=NaN;
            zeta_u=0.5*(zeta(1:len_lat_rho-1,:)+zeta(2:len_lat_rho,:));
            zeta_v=0.5*(zeta(:,1:len_lon_rho-1)+zeta(:,2:len_lon_rho));
            g.depth=g.h+zeta;
            g.depth_u=g.h_u + zeta_u;
            g.depth_v=g.h_v + zeta_v;

    %         s_max = 34.573;


    %         % size of grids
    %         mask3d_rho=repmat(g.mask_rho,[1 1 g.N]);
    %         mask3d_rho=permute(mask3d_rho,[3 1 2]);
    % 
    %         [r1,c1] = size ( g.lon_u );
    %         mask3d_u=repmat(g.mask_u,[1 1 g.N]);
    %         mask3d_u=permute(mask3d_u,[3 1 2]);
    % 
    %         [r2,c2] = size ( g.lon_v );
    %         mask3d_v=repmat(g.mask_v,[1 1 g.N]);
    %         mask3d_v=permute(mask3d_v,[3 1 2]);

            % transport from surface to which depth (m)
            selectdepth=5000;
            if ( selectdepth > 0 )
                selectdepth = selectdepth*-1;
            end
            maxdepth = 5000;

            sn=0;
            disp(['month ',': ',num2str(month)])

            hold on

            for st=1:t_point
    %             disp([char(point_name(st)),' Depth to ',num2str(maxdepth)])
                endpt_lon(1) = kts(sn+1);
                endpt_lat(1) = kts(sn+2);
                endpt_lon(2) = kts(sn+3);
                endpt_lat(2) = kts(sn+4);
                sn=sn+4;

    % %             [indw, inde, inds, indn] = [corner_endpt_row(1), corner_endpt_row(2), corner_endpt_col(1), corner_endpt_col(2)]
    %             [corner_endpt_col(1), corner_endpt_col(2), corner_endpt_row(1), corner_endpt_row(2)] = ...
    %                 findind_Y(0,[endpt_lon(1),endpt_lon(2),endpt_lat(1),endpt_lat(2)],g.lon_rho', g.lat_rho');

                for cpts=1:2 % corner points
            % %         find points that are close with stations 
                    dist = sqrt(  ( g.lon_rho - endpt_lon(cpts) ).*( g.lon_rho - endpt_lon(cpts) ) + ...
                        ( g.lat_rho - endpt_lat(cpts) ).*( g.lat_rho - endpt_lat(cpts) ) );
                    ind=find( min( dist(:) ) == dist );
                    % closest points row and column indice
                    row_index = mod ( ind - 1, len_lat_rho ) + 1;
                    col_index = floor( (ind - 1) / len_lat_rho ) + 1;
                    corner_endpt_col(cpts)=col_index(1);
                    corner_endpt_row(cpts)=row_index(1);
                end


                % my xy_transport_onefile works only if corner_endpt_row(2) < corner_endpt_row(1).
                if( corner_endpt_row(2) > corner_endpt_row(1)  )
                    tmp_col=corner_endpt_col(2);
                    tmp_row=corner_endpt_row(2);
                    corner_endpt_col(2)=corner_endpt_col(1);
                    corner_endpt_row(2)=corner_endpt_row(1);
                    corner_endpt_col(1)=tmp_col;
                    corner_endpt_row(1)=tmp_row;
    %                 beep
                    disp(' === switching two end points === ')
                end


                % longitude and latitude coordinate.

                for i=1:length(endpt_lat)
                    xx(i)=g.lon_rho(corner_endpt_row(i),corner_endpt_col(i));
                    yy(i)=g.lat_rho(corner_endpt_row(i),corner_endpt_col(i));
                end
                distance_r = m_lldist ( xx, yy )*1000;

               %  transect information

                if( corner_endpt_col(2) >= corner_endpt_col(1) )
                    delj=1;
                    west2east_transect=1; % previously zonaltransect
                else
                    delj=-1;
                    west2east_transect=0; % previously meridionaltransect
                end

                if( corner_endpt_row(2) > corner_endpt_row(1) )
                    deli=1;
                else
                    deli=-1;
                end

                xzero=g.lon_rho( corner_endpt_row(1), corner_endpt_col(1) );
                yzero=g.lat_rho( corner_endpt_row(1), corner_endpt_col(1) );
                xone=g.lon_rho( corner_endpt_row(2), corner_endpt_col(2) );
                yone=g.lat_rho( corner_endpt_row(2), corner_endpt_col(2) );
                slope=( yone-yzero) / (xone - xzero);
                % A x0 + B y0 + C = 0;
                A=slope;
                B=-1;
                C=-slope*xzero+yzero;
                D=sqrt( A*A + B*B );
                % distance between point and line = abs( A x + B y + C ) / D. 

                %   grid information

                avgdxdy = mean([ mean( mean( g.dx ) ), mean( mean( g.dy ) ) ]);
                maxdxdy = max([ max( max( g.dx ) ), max( max( g.dy ) ) ]);

                ynum = 0; mnum=0;
                ynum=ynum+1;
                mnum=mnum+1;

                ubar=ncread(filename, 'ubar')';
                ubar(ubar>=1000)=NaN;
                ubar(ubar<=-1000)=NaN;
                vbar=ncread(filename, 'vbar')';
                vbar(vbar>=1000)=NaN;
                vbar(vbar<=-1000)=NaN;
    % 
    %             icount=1;
    %             col_index(icount)=corner_endpt_col(1); % yind of the present point
    %             row_index(icount)=corner_endpt_row(1); % xind of the present point
    %             vpoint=1; % flag of the vpoint u==0 or v==1)
    %             if (deli==0)
    %                 vpoint=0;
    %             end
    %             on_vpoint(icount)=vpoint; % flag of the point (u==0 or v==1), accumulated
    %             xpoint=g.lon_u( row_index(icount), col_index(icount) );  % lon of the present point
    %             ypoint=g.lat_u( row_index(icount), col_index(icount) );  % lat of the present point
    % 
    %             signline(icount)=(ypoint - yzero) - slope*(xpoint -xzero); % sign of the line, if partial slope is bigger than the slope, sign is positive.
    %             dist(icount)= abs( A*xpoint + B*ypoint + C ) / D;  % distance from slope line
    %             tmp_dist=dist(icount);
    %             dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000; % distance between start and end point
    %             flag_approach = 1; % if flag_approach is 1, continue to calculate

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
                    dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                    flag_approach = 1;

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
    %                                 plot(  xpoint,  ypoint , 'ro')
                                    dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                                else
                                    tmp_dist=dist(icount);
                                    vpoint=0;
                                    icount=icount-1;
    %                                 plot(  xpoint,  ypoint , 'gx')
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
    %                                 plot(  xpoint,  ypoint , 'ro')
                                    dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                                else
                                    tmp_dist=dist(icount);
                                    vpoint=1;
                                    icount=icount-1;
    %                                 plot(  xpoint,  ypoint , 'gx')
                                end

                            end % if ( on_vpoint == 1 )

                            if( icount > 3 &&  dist2endpoint(icount) > dist2endpoint(icount-3) )
                                flag_approach = 0;
                            end

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
    %                                 plot(  xpoint,  ypoint , 'ro')
                                    dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                                else
                                    tmp_dist=dist(icount);
                                    vpoint=0;
                                    icount=icount-1;
    %                                 plot(  xpoint,  ypoint , 'gx')
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
    %                                 plot(  xpoint,  ypoint , 'ro')
                                    dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                                else
                                    tmp_dist=dist(icount);
                                    vpoint=1;
                                    icount=icount-1;
    %                                 plot(  xpoint,  ypoint , 'gx')
                                end

                            end % if ( on_vpoint == 1 )

                            if( icount > 3 &&  dist2endpoint(icount) > dist2endpoint(icount-3) )
                                flag_approach = 0;
                            end
                        end % while

                    end % if (west2east_transect)

    %             total_temp=0;
                total_volume=0;

                for index=1:icount

                    vpoint=on_vpoint(index);

    % % %                 calculate transport using xy_transport_function
                    if ( vpoint ==  1 )  
                        jj=row_index(index);
                        i1=col_index(index);
                        i2=col_index(index);
                        % actual calculation of transport
                        sum_segment=0;
    %                     end
                        for i=i1:i2
                            temp_segment= g.depth_v(jj,i)*vbar(jj,i)*g.dx_v(jj,i);
                            if isfinite(temp_segment)
                                sum_segment=sum_segment+g.depth_v(jj,i)*vbar(jj,i)*g.dx_v(jj,i);
                            end
                        end

                    elseif  ( vpoint == 0 )  
                        ii=col_index(index); 
                        j1=row_index(index);
                        j2=row_index(index);
                        % actual calculation of transport
                        sum_segment=0;
                        for j=j1:j2
                            temp_segment= g.depth_u(j,ii)*ubar(j,ii)*g.dy_u(j,ii);
                            if isfinite(temp_segment)
                                sum_segment=sum_segment+g.depth_u(j,ii)*ubar(j,ii)*g.dy_u(j,ii);
                            end
                        end
                    end

                    if (isfinite(sum_segment)==1)
                        if( west2east_transect == 0 && vpoint == 1 )
                            total_volume = total_volume - sum_segment;
                        else
                            total_volume = total_volume + sum_segment;
                        end
                    end
               end
               disp(['transport ',char(point_name(st)),': ',num2str(total_volume/1e+6)])
               trans(mnum,st)    = total_volume/1e+6; 
            end  
            hold off
            fid = fopen(outfile,'a+');
        %     fprintf(fid,'%8.3f %8.3f %8.3f  \n', trans');
            fprintf(fid,'%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n', trans');
            fclose(fid);
        end
    end

end
