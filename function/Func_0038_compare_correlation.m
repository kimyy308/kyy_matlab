function p = Func_0038_compare_correlation(r1,r2,n1,n2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% compare correlation
% ref: Cohen, J., P. Cohen, S. G. West, and L. S. Aiken, 2013 
% Applied multiple regression/correlation analysis for the behavioral sciences.  Routledge.
% 49 ~ 50 pages
%
%  input:
%  r1         correlation 1
%  r2         correlation 2
%  n1         # of samples for correlation 1
%  n2         # of samples for correlation 2
%
%  output:
%  p          p value, the probability that H0 
%            (the correlation coefficiets are not different) is correct
%
%  e-mail:      kimyy308@pusan.ac.kr
%
%  Based on https://kr.mathworks.com/matlabcentral/fileexchange/44658-compare_correlation_coefficients
%
%  Updated      07-Nov-2024 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_r1 = 0.5*log((1+r1)./(1-r1));
t_r2 = 0.5*log((1+r2)./(1-r2));
z = (t_r1-t_r2)./sqrt(1./(n1-3)+1./(n2-3));
p = (1-normcdf(abs(z),0,1))*2;


end
