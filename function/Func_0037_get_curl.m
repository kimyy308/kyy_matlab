function [curlZ] = Func_0037_get_curl(lon,lat,Tx,Ty)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [lv, pc, var_exp] = Func_0037_get_curl(lon,lat,X);
%
% calculate curl from 2-d or 3-d field [x,y,t]
%
%  input:
%  lon         longitude [x,y]
%  lat         latitude  [x,y]
%  U           zonal velocity       [x,y] or [x,y,t]
%  V           meridional velocity  [x,y] or [x,y,t]
%
%  output:
%  curlx       curl or vorticity
%
%  e-mail:      kimyy308@pusan.ac.kr
%
%  Based on Ramkrushn S. Patel's code (ra_windstrcurl.m; ramkrushn.scrv89@gmail.com)
%  must not use curvlinear grid. normal roms grid is assumed
%  Updated      31-May-2024 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % computation of curl
% [lt, ln]=size(u);
% a=diff(lat); a = fix(a*10^4)/10^4;
% aa=NaN*ones(length(a)-1,1);
% for ii=1:length(a)-1
%     if (a(ii) == a(ii+1))
%         aa(ii)=a(ii);
%     else
%         error('Latitude difference is not consistance')
%     end % endif
%     dlat=mean(aa);
% end % endfor
% clear ii
% deltay=dlat*111176;
% curlZ=NaN(lt, ln);
% long=NaN(lt, ln);
% for ii=1:lt
%     for jj=1:ln
%         long(ii, jj)=lon(jj)*111176*cos(lat(ii)*rad);
%         %long(i,j)=lon(j)*6378137*rad*cos(lat(i)*rad);
%         % [m] earth radious in meters= 6,378,137.0 m.. from wikipedia.
%     end % endfor
% end % endfor
% clear ii jj

Tx(Tx>10e6)=NaN;
Ty(Ty>10e6)=NaN;
Tx(Tx<-10e6)=NaN;
Ty(Ty<-10e6)=NaN;

if ndims(Tx)==2

    [ln, lt]=size(Tx);
    
    dlat(1,:)=m_lldist(repmat(lon(1,1),[1,length(lat(1,:))]) ,lat(1,:),length(lat(1,:)));
    dy=repmat(dlat, [size(lon,1),1]); % unit: kilometer
    dy(:,end+1)=dy(:,end);
    
    for li=1:length(lon(1,:))
        dlon(1,li)=m_lldist([lon(1,li), lon(2,li)], [lat(1,li), lat(2,li)]);
    end
    dx=repmat(dlon, [size(lon,1),1]); % unit: kilometer
    
    
    
    % Centeral difference method in x and y
    for ii=2:ln-1
        for jj=2:lt-1
            curlZ(ii, jj)=(Ty(ii, jj+1)-Ty(ii, jj-1))/((dx(ii, jj+1)+dx(ii, jj-1))) - ...
                (Tx(ii+1, jj)-Tx(ii-1, jj))/(2*dy(ii,jj)) ;
        end % endfor
    end % endfor
    clear ii jj
    % Forward difference method in x and y
    for jj=1:lt-1
        curlZ(1, jj)=(Ty(1, jj+1)-Ty(1, jj))/(dx(1, jj+1)+dx(1, jj))/2 - ...
            (Tx(2, jj)-Tx(1, jj))/dy(1, jj) ;
    end
    for ii=1:ln-1
        curlZ(ii, 1)=(Ty(ii, 2)-Ty(ii, 1))/(dx(ii, 2)+dx(ii, 1))/2 - ...
            (Tx(ii, 2)-Tx(ii, 1))/dy(ii, 1) ;
    end
    clear ii jj
    curlZ(1, ln)=curlZ(1, ln-1);
    % Backward difference method in x and y
    for ii=2:ln
        curlZ(ii, lt)=(Ty(ii, lt)-Ty(ii, lt-1))/(dx(ii, lt)+dx(ii, lt-1))/2 - ...
            (Tx(ii, lt)-Tx(ii-1, lt))/dy(ii,lt) ;
    end
    for jj=2:lt-1
        curlZ(ln, jj)=(Ty(ln, jj)-Ty(ln, jj-1))/(dx(ln, jj)+dx(ln, jj-1))/2 - ...
            (Tx(ln, jj)-Tx(ln-1, jj))/dy(ln,jj) ;
    end
    clear ii jj
    % curlZ(lt, 1)=curlZ(lt, lt-1);


else

    for ti=1:size(Tx,3)
        [ln, lt, ltime]=size(Tx);
        
        dlat(1,:)=m_lldist(repmat(lon(1,1),[1,length(lat(1,:))]) ,lat(1,:),length(lat(1,:)));
        dy=repmat(dlat, [size(lon,1),1]); % unit: kilometer
        dy(:,end+1)=dy(:,end);
        
        for li=1:length(lon(1,:))
            dlon(1,li)=m_lldist([lon(1,li), lon(2,li)], [lat(1,li), lat(2,li)]);
        end
        dx=repmat(dlon, [size(lon,1),1]); % unit: kilometer
        
        
        
        % Centeral difference method in x and y
        for ii=2:ln-1
            for jj=2:lt-1
                curlZ(ii, jj, ti)=(Ty(ii, jj+1, ti)-Ty(ii, jj-1, ti))/((dx(ii, jj+1)+dx(ii, jj-1))) - ...
                    (Tx(ii+1, jj, ti)-Tx(ii-1, jj, ti))/(2*dy(ii,jj)) ;
            end % endfor
        end % endfor
        clear ii jj
        % Forward difference method in x and y
        for jj=1:lt-1
            curlZ(1, jj, ti)=(Ty(1, jj+1, ti)-Ty(1, jj, ti))/(dx(1, jj+1)+dx(1, jj))/2 - ...
                (Tx(2, jj, ti)-Tx(1, jj, ti))/dy(1, jj) ;
        end
        for ii=1:ln-1
            curlZ(ii, 1, ti)=(Ty(ii, 2, ti)-Ty(ii, 1, ti))/(dx(ii, 2)+dx(ii, 1))/2 - ...
                (Tx(ii, 2, ti)-Tx(ii, 1, ti))/dy(ii, 1) ;
        end
        clear ii jj
        curlZ(1, ln)=curlZ(1, ln-1);
        % Backward difference method in x and y
        for ii=2:ln
            curlZ(ii, lt, ti)=(Ty(ii, lt, ti)-Ty(ii, lt-1, ti))/(dx(ii, lt)+dx(ii, lt-1))/2 - ...
                (Tx(ii, lt, ti)-Tx(ii-1, lt, ti))/dy(ii,lt) ;
        end
        for jj=2:lt-1
            curlZ(ln, jj, ti)=(Ty(ln, jj, ti)-Ty(ln, jj-1, ti))/(dx(ln, jj)+dx(ln, jj-1))/2 - ...
                (Tx(ln, jj, ti)-Tx(ln-1, jj, ti))/dy(ln,jj) ;
        end
        clear ii jj
        % curlZ(lt, 1)=curlZ(lt, lt-1);
    end

end


end
