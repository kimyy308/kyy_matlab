clear all;close all;clc

for year=2009:2010
    for month=1:1:12

        nc=netcdf(['SODA_2.2.4_',num2str(year),num2char(month,2),'.cdf'],'write');
        SSH=nc{'SSH'}(:);
        SALT=nc{'SALT'}(:);
        TEMP=nc{'TEMP'}(:);
        U=nc{'U'}(:);
        V=nc{'V'}(:);
        [z,m,n]=size(U);
        for k=1:1:z
            for i=1:1:m
                for j=1:1:n
                    if SSH(i,j) < -999
%                         SALT(k,i,j)=NaN;
%                         TEMP(k,i,j)=NaN;
%                         U(k,i,j)=NaN;
%                         V(k,i,j)=NaN;
                        SSH(i,j)=NaN;
                    end
                end
            end
        end
        nc{'SSH'}(:)=SSH;
%         nc{'SALT'}(:)=SALT;
%         nc{'TEMP'}(:)=TEMP;
%         nc{'U'}(:)=U;
%         nc{'V'}(:)=V;
        close(nc)
    end
end
