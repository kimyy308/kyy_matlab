function coef = velocity_correlation(u_rho_1,v_rho_1,u_rho_2,v_rho_2)


    u1 = u_rho_1(logical(~isnan(u_rho_1).*~isnan(u_rho_2).*~isnan(v_rho_1).*~isnan(v_rho_2)));
    v1 = v_rho_1(logical(~isnan(u_rho_1).*~isnan(u_rho_2).*~isnan(v_rho_1).*~isnan(v_rho_2)));
    u2 = u_rho_2(logical(~isnan(u_rho_1).*~isnan(u_rho_2).*~isnan(v_rho_1).*~isnan(v_rho_2)));
    v2 = v_rho_2(logical(~isnan(u_rho_1).*~isnan(u_rho_2).*~isnan(v_rho_1).*~isnan(v_rho_2)));

%     u1 = [ 1 3 2 4 1 3 2 4];
%     v1 = -[ 1 2 3 4 5 6 7 8];
%     u2 = [ 1 4 3 4 1 1 3 4];
%     v2 = -[ 1 2 2 5 1 2 2 1];

    uu_cov=cov(u1,u2);
    vv_cov=cov(v1,v2);
    uv_cov=cov(u1,v2);
    vu_cov=cov(v1,u2);

% % % method 1
%     
%     sigma_1=double(inv(cov(u1,v1)));
% 
%     sigma_2(1,1)=uu_cov(1,2);
%     sigma_2(2,2)=vv_cov(1,2);
%     sigma_2(2,1)=uv_cov(2,1);
%     sigma_2(1,2)=uv_cov(2,1);
%     
%     sigma_3=double(inv(cov(u2,v2)));
% 
%     sigma_4(2,1)=uv_cov(1,2);
%     sigma_4(1,2)=uv_cov(1,2);
%     sigma_4(1,1)=uu_cov(2,1);
%     sigma_4(2,2)=vv_cov(2,1);

%     trace(sigma_1 * sigma_2 * sigma_3 * sigma_4)

% % method 2
    cov_sig_1=cov(u1,v1);
    u1u1=cov_sig_1(1,1);
    u1v1=cov_sig_1(1,2);
    v1u1=cov_sig_1(2,1);
    v1v1=cov_sig_1(2,2);
    cov_sig_2=cov(u2,v2);
    u2u2=cov_sig_2(1,1);
    u2v2=cov_sig_2(1,2);
    v2u2=cov_sig_2(2,1);
    v2v2=cov_sig_2(2,2);
    u1u2=uu_cov(1,2);
    u2u1=uu_cov(2,1);
    v1v2=vv_cov(1,2);
    v2v1=vv_cov(2,1);
    u1v2=uv_cov(2,1);
    v2u1=uv_cov(1,2);
    u2v1=vu_cov(2,1);
    v1u2=vu_cov(1,2);

    f= u1u1 * (u2u2*v1v2^2 + v2v2*v1u2^2) ...
        + v1v1 * (u2u2*u1v2^2 + v2v2*u1u2^2) ...
        + 2*(u1v1*u1v2*v1u2*u2v2) ...
        + 2*(u1v1*u1u2*v1v2*u2v2) ...
        - 2*(u1u1*v1u2*v1v2*u2v2) ...
        - 2*(v1v1*u1u2*u1v2*u2v2) ...
        - 2*(u2u2*u1v1*u1v2*v1v2) ...
        - 2*(v2v2*u1v1*u1u2*v1u2);
    g = (u1u1*v1v1 - u1v1^2) * (u2u2*v2v2-u2v2^2);
    coef=sqrt(f/g/2.0)
    
    
% % % method 3
%     u1u1 = cov(u1,u1);
%     u1u2 = cov(u1,u2);
%     u2u1 = cov(u2,u1);
%     u2u2 = cov(u2,u2);
%     v1v1 = cov(v1,v1);
%     v1v2 = cov(v1,v2);
%     v2v1 = cov(v2,v1);
%     v2v2 = cov(v2,v2);
%     u1v1 = cov(u1,v1);
%     v1u1 = cov(v1,u1);
%     u1v2 = cov(u1,v2);
%     v2u1 = cov(v2,u1);
%     u2v1 = cov(u2,v1);
%     v1u2 = cov(v1,u2);
%     u2v2 = cov(u2,v2);
%     v2u2 = cov(v2,u2);
%     
%     f= u1u1 * (u2u2*v1v2^2 + v2v2*v1u2^2) ...
%         + v1v1 * (u2u2*u1v2^2 + v2v2*u1u2^2) ...
%         + 2*(u1v1*u1v2*v1u2*u2v2) ...
%         + 2*(u1v1*u1u2*v1v2*u2v2) ...
%         - 2*(u1u1*v1u2*v1v2*u2v2) ...
%         - 2*(v1v1*u1u2*u1v2*u2v2) ...
%         - 2*(u2u2*u1v1*u1v2*v1v2) ...
%         - 2*(v2v2*u1v1*u1u2*v1u2);
%     g = (u1u1*v1v1 - u1v1^2) * (u2u2*v2v2-u2v2^2);
%     coef=sqrt(f/g/2.0)
    
% %     perfect correlation example
%  1. same vector series
%  2. constant value is multiplied
%  3. rotated by a constant angle
%  4. both multiplied by a constant and rotated by a constant angle
return