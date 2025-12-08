function [L, U, p, Rflag] = Func_0039_compare_correlation_Siegert(r1,r2, r12, n1,n2, alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% compare correlation
% ref: Siegert et al., 2017 
%
%  input:
%  r1         correlation 1
%  r2         correlation 2
%  n1         # of samples for correlation 1
%  n2         # of samples for correlation 2
%
%  output:
%  L          lower boundary of insignificant value
%  U          upper boundary of insignificant value
%  p          p-value by one-sided T-test based on T2 from Siegert et al., 2017
%  Rflag      When R is negative (undefined situation, but here p replaced with insignificant value 1).
%            
%
%  e-mail:      kimyy308@pusan.ac.kr
%
%
% it is based on null hypothesis (correlation difference r2-r1 = 0),
% so D > U, D<L means significantly different.
%
%  Updated      18-Apr-2025 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% upper & lower significantly consistent value
    z_ref_lower=norminv(alpha/2);
    z_ref_upper=norminv(1-alpha/2);
    
    z1=0.5 * log((1+r1)./(1-r1));
    z2=0.5 * log((1+r2)./(1-r2));
    
    l1=tanh(z1 + z_ref_lower/sqrt(n1-3));
    u1=tanh(z1 + z_ref_upper/sqrt(n1-3));
    l2=tanh(z2 + z_ref_lower/sqrt(n2-3));
    u2=tanh(z2 + z_ref_upper/sqrt(n2-3));
    
    c_ab = ((r12-0.5*r1*r2)*(1-r1^2-r2^2-r12^2)+r12^3) / ((1-r1^2)*(1-r2^2));
    
    L= (r2-r1) - sqrt( (r2-l2)^2 + (u1-r1)^2 - 2*c_ab*(r2-l2)*(u1-r1) );
    U= (r2-r1) + sqrt( (u2-r2)^2 + (r1-l1)^2 - 2*c_ab*(u2-r2)*(r1-l1) );
    
    %% test statistics (T test for overlapped correlation difference (r_ay, r_by))
    R = (1 - r1^2 - r2^2 - r12^2) + (2*r1*r2*r12);
    
    % p = 2 * (1 - tcdf(abs(T), n1-3)); %two_sided test
    % one sided test
    
    Rflag = 0;
    if R>=0
        T= (r2-r1) * sqrt( ((n1-1)*(1+r12)) / (2*(n1-1)/(n1-3)*R + 1/4*(r1+r2)^2*(1-r12)^3));
        if (r2-r1)>=0
            p = 1 - tcdf(T, n1-3);
        else
            p  = tcdf(T, n1-3);
        end
    elseif isnan(R)==1
        p = NaN;
        Rflag = NaN;
    % else % normal fisher test
    %     disp('caution, R<0 occured')
    %     t_r1 = 0.5*log((1+r1)./(1-r1));
    %     t_r2 = 0.5*log((1+r2)./(1-r2));
    %     z = (t_r1-t_r2)./sqrt(1./(n1-3)+1./(n2-3));
    %     p = (1-normcdf(abs(z),0,1))*2;
    else % undefined due to R < 0
    %     disp(['caution, R<0 occured, ', 'R=',num2str(R), ', r1=', num2str(r1), ', r2=', num2str(r2), ', r12=', num2str(r12)]);
        p = 1;
        Rflag = 1;
    end

end
