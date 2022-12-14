function [lv, pc, var_exp] = Func_0024_EOF_3d(data,X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [lv, pc, var_exp] = Func_0024_EOF_3d(data,X);
%
% Empirical Orthogonal Function (EOF) using 3-d data[x,y,t]
%
%  input:
%  X            The number of modes that you'd like to get (integer)
%  data         The 3-dimensional input data [x,y,t] to be analyzed
%
%  output:
%  lv           Dimensionless Eigenvector(Loading vector) [x,y,l]
%  pc           Eigenvalue(Principal Component (PC)) with dimension [t,l]
%  var_exp      Explained variance(%) [l]
%
%  e-mail:      kimyy308@pusan.ac.kr
%
%  Updated      14-Dec-2022 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% get the original size of data
    [M,N,L]=size(data);
    
    %% get grids of finite value
    mask_finite=isfinite(sum(data,3)); % get finite mask considering whole period
    K=sum(mask_finite(:)); % length of valid values at certain time
    
    %% Create 2-D matrix by finite values only
    data_compact=zeros(K,L);
    for i=1:L
        tmpdata=data(:,:,i);
        data_compact(:,i)=tmpdata(mask_finite);
    end
    
    %% Singular Value Decomposition
    [u,s,pc]= svd(data_compact,0);
    
    %% get loading vector
    u=u(:,1:X); % get valid modes
    lv=NaN(M,N,X); % lv initialization
    for i=1:X
        tmpdata=NaN(M,N);
        tmpdata(mask_finite)=u(:,i); % recast u to the raw grid
        lv(1:M,1:N,i)=tmpdata;
    end
    
    %% get principal component(pc)
    pc=pc(:,1:X); % get valid modes
    
    %% get explained variance (%)
    var=diag(s).^2;
    var_exp=var/sum(var) .*100; % get percentage
    var_exp=var_exp(1:X); % get valid modes

end


% [M,N,L]=size(SSH);
% SSH_2d=reshape(SSH,[M*N, L]);
% mask_nan=isnan(SSH_2d);
% SSH_2d(mask_nan)=0;
% [u,s,pc]= svd(SSH_2d,0);
% var=diag(s).^2;
% var_exp=var/sum(var);
% plot(var_exp)
% u(mask_nan)=NaN;
% lv=reshape(u, [M N L]);
% pcolor(lv(:,:,1)'.*pc(1,1)); shading flat; colorbar;



% % % The EOF program
% % %
% % % [seof,teof, pecenv]=EOF_PAN(m,n,l,inputmat)
% % % 
% % % seof(m,n,l) is Spatial EOF components,  where l is EOF mode number, 
% % %
% % % m and n are the spatial coordinates
% % %
% % % teof(l,X) is  temporal EOF components, where X is mode number
% % %
% % % df(m,n,l) is a smaple matrix, where m and are spatial coordinates
% % %
% % % l is time series number
% % %
% % % pecenv is pecentage of the contribution of every EOF component to the total variance.
% % 
% % % Copyright reserved, Jan., 2000
% % 
% % 
% % 
% % 
% % function  [seof,teof, pecenv, s]=EOF_PAN(m,n,l,inputmat)
% % 
% % ln(1:m,1:n)=0;
% % 
% % for i=1:l
% % ln=ln+inputmat(:,:,i);
% % end
% % 
% % ln(not(isnan(ln)))=0;
% % 
% % 
% % k=0;
% % 
% % for i=1:m
% %     for j=1:n
% %     
% %     if not(isnan(ln(i,j)))
% %     
% %     k=k+1;
% %     z(1:l)=inputmat(i,j,:);
% %     samplemat1(1:l,k)=z';
% %     clear z
% %     
% %     end
% %     
% %     end
% % end
% % 
% % 
% % samplemat2=samplemat1-ones(l,1)*mean(samplemat1);
% % 
% % workmat=samplemat2';    
% % 
% % 
% %     
% %    [u,s,v]=svd(workmat,0);
% %    
% %      
% %   k1=0;
% % for i=1:m
% %    for j=1:n
% %     
% %      if not(isnan(ln(i,j)))
% %       k1=k1+1;
% %         seof(i,j,1:l)=u(k1,:);
% %         
% %        else
% %        
% %         seof(i,j,1:l)=NaN; 
% %         
% %         
% %      end
% %    
% %      end
% %  end   
% %   
% %  teof=v;
% %  
% %  variance1=diag(s);
% %  
% %  variance=variance1.*variance1;
% %  
% %  pecenv=variance/sum(variance);
