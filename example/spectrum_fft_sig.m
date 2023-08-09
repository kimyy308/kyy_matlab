clc; close all; clear all;
% filename='/Volumes/kyy_raid/kimyy/Model/CMIP/CMIP6_analysis/tos/NA_reg/feb/runmean_9/runm_tos_NA_reg_ACCESS-CM2_feb.nc';
filename='/Volumes/kyy_raid/kimyy/Model/CMIP/CMIP6_analysis/tos/NA_reg/feb/runmean_9/runm_tos_NA_reg_MIROC6_feb.nc';

ncdisp(filename);
tos=squeeze(ncread(filename, 'tos'));
tos=tos(8:57);
n=length(tos);
nyquist=1/2;
freq=(1:n/2)/(n/2)*nyquist;
period=1./freq;
pft=fft(tos);
pft=pft;
pft(1) = [];
% pft_real=real(pft);
% pft_imag=imag(pft);
pxx=abs(pft(1:floor(n/2))).^2 ./n .*2;
% pxx=abs(pft(1:floor(n/2)));
% pxx=pxx./n;
% pxx=pxx.*2.0;
% for i=1:length(pxx)
%     f(i)=1/(N/i);
% end
r=autocorr(tos,2);
r1=r(2);
% Chi=5.991; %Chi-Square values 95% (chi2inv(0.95, 2))
Chi=chi2inv(0.95,2);
for i=1:floor(n/2)
    Red(i)=(1-r1^2)/(1+r1^2-2*r1*cos(2*pi*freq(i)))/2*Chi; % Rednoise significance level
end
plot(freq,pxx)
hold on
plot(freq,Red)
hold off
% xlim([0 0.1])
ylim([0 70])

[psd, psd_sig, freq_2, var_2]=Func_0027_get_PSD_siglev(tos);

figure;
plot(freq,psd)
hold on
plot(freq,psd_sig)
hold off
% xlim([0 0.1])
ylim([0 70])

% plot(1./freq,pxx,'-*')
% hold on
% plot(1./freq,Red,'-o')
% hold off
% xlim([9 30])
% % % % % 
% % % % % % real(pft)
% % % % % % imag(pft)
% % % % % 
% % % % % % plot(period,pxx)
% % % % % % xlim([2 30])
% % % % % 
% % % % % 
% % % % % % % %         fu2=fft(tg_v2);
% % % % % % % %         fu1(1) = [];
% % % % % % % %         fu2(1) = [];
% % % % % % % %         power1=abs(fu1(1:floor(n/2))).^2;
% % % % % 
% % % % % 
% % % % % % plot(tos);
% % % % % 
% % % % % % [pxx,w]=periodogram(tos(8:57), 'power', 'ConfidenceLevel', 0.95);
% % % % % % [pxx,w]=periodogram(tos, 'power');
% % % % % % per=1.0./w*50.0/pi;
% % % % % % plot(per,pxx)
% % % % % % xlim([15 27])
% % % % % 
% % % % % % fft(tos)
% % % % % % 1961/2010 (8~57)
% % % % % % 2047/2096 (94~143)
% % % % % % 1954:2096
% % % % % 
% % % % % % [pxx,w]=pspectrum(tos(8:57))
% % % % % 
% % % % % % fftshift(fft(tos))/143
% % % % % 
% % % % % % [f_final,s_final,phase]=fftspectrum(tos,1:143);
% % % % % 
% % % % % 
% % % % % 
% % % % % %% code from hk lee
% % % % % Data=tos;
% % % % % % Data=pxx(2:end);
% % % % % 
% % % % % N=length(Data);
% % % % % Y_bar=0;
% % % % % Y_std=0;
% % % % % for i= 1:N
% % % % %     Y_bar=Y_bar+Data(i);
% % % % % end
% % % % % 
% % % % % 
% % % % %  Y_bar=Y_bar/N;
% % % % %  for i=1 :N
% % % % %      Y_std=Y_std+(Data(i)-Y_bar)^2;
% % % % %  end
% % % % %  Y_std=sqrt(Y_std/(N-1));
% % % % % 
% % % % % for k=1:N/2
% % % % %     A(k)=0;
% % % % %     B(k)=0;
% % % % %    
% % % % % for t=1: N
% % % % %    
% % % % %     A(k)=A(k)+((2/N)*Data(t)*cos(2*pi*k*t/N));
% % % % %     B(k)=B(k)+((2/N)*Data(t)*sin(2*pi*k*t/N));
% % % % %    
% % % % % end
% % % % % 
% % % % % 
% % % % % C(k)=sqrt((A(k)^2)+(B(k)^2));
% % % % % 
% % % % % 
% % % % % end
% % % % % for i=1:length(C)
% % % % %     f(i)=1/(N/i);
% % % % % end
% % % % % r=autocorr(Data,2);
% % % % % r1=r(2);
% % % % % Chi=5.991; %Chi-Square values 95%
% % % % % for i=1:N/2
% % % % % Red(i)=(1-r1^2)/(1+r1^2-2*r1*cos(2*pi*f(i)))/2*Chi; % Rednoise significance level
% % % % % end
% % % % % 
% % % % % % 이 매트랩 코드쓰시면 됩니다. 첨부파일에 3번식이 Rednoise significance level 입니다.
% % % % % % Data에 아무 데이터 값 넣으시면 되고
% % % % % % plot할때
% % % % % 
% % % % % figure
% % % % % plot(f,C.^2*N)
% % % % % hold on
% % % % % plot(f,Red)
% % % % % hold off
% % % % % % 로 비교하시면됩니다
% % % % % % 그리고 1분 데이터일 때 f/60
% % % % % % 10초 데이터일 때 f/10 입니다.
% % % % % % 
% % % % % % 
% % % % % 
% % % % % 
% % % % % % % load sunspot.dat
% % % % % % % % relNums=sunspot(:,2);
% % % % % % % relNums=tos;
% % % % % % % [pxx,f] = periodogram(relNums,[],[],1);
% % % % % % % plot(1/f,10*log10(pxx))
% % % % % % % % xlabel('Cycles/Year')
% % % % % % % xlabel('Year')
% % % % % % % ylabel('dB / (Cycles/Year)')
% % % % % % % title('Periodogram of Relative Sunspot Number Data')
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % % stlee's code
% % % % % 
% % % % % % % % clear all;close all;clc
% % % % % % % % 
% % % % % % % % exp_flag = 1; % 1:NODA / 2:EnOI / 3:EnKF
% % % % % % % % dim_flag = 2; % 1: point / 2: enitre domain
% % % % % % % % regional_flag = 1;
% % % % % % % % 
% % % % % % % % switch exp_flag
% % % % % % % %     case 1
% % % % % % % %         var_name1='NODA';
% % % % % % % %         var_name2='noda';
% % % % % % % %     case 2
% % % % % % % %         var_name1='EnOI';
% % % % % % % %         var_name2='enoi';
% % % % % % % %     case 3
% % % % % % % %         var_name1='EnKF';
% % % % % % % %         var_name2='enkf';
% % % % % % % % end
% % % % % % % % 
% % % % % % % % grdpath = 'D:\model_output\DA\NODA\2015\';
% % % % % % % % grdname = [grdpath,'ocean_his_0001.nc'];
% % % % % % % % 
% % % % % % % % lon = ncread(grdname,'lon_rho');
% % % % % % % % lat = ncread(grdname,'lat_rho');
% % % % % % % % rmask = ncread(grdname,'mask_rho');
% % % % % % % % rmask(find(rmask==0))=nan;
% % % % % % % % 
% % % % % % % % if (exp_flag ==2)
% % % % % % % %     filename1 = 'D:\model_output\DA\EnOI\enoi_ensemble_mean.mat';
% % % % % % % %     d1=load(filename1);
% % % % % % % %     v1=d1.v1_m; % sst
% % % % % % % %     v2=d1.v2_m; % ssha
% % % % % % % % else
% % % % % % % %     filename1 = strcat('D:\model_output\DA\',var_name1,'\2015\',var_name2,'_sst.mat');
% % % % % % % %     filename2= strcat('D:\model_output\DA\',var_name1,'\2015\',var_name2,'_zeta.mat');
% % % % % % % % 
% % % % % % % %     d1 = load(filename1);
% % % % % % % %     d2 = load(filename2);
% % % % % % % %     v1 = d1.sst; % sst(noda)
% % % % % % % %     v2 = d2.zeta; % ssha(noda)
% % % % % % % % end
% % % % % % % % for i=1:365
% % % % % % % %     v2(:,:,i) = v2(:,:,i) - nanmean(v2(:,:,i),[1 2]);
% % % % % % % % end
% % % % % % % % 
% % % % % % % % t1 = datenum(2018,1,1):datenum(2018,12,31);
% % % % % % % % t_tick = datenum(2018,[1:12],1);
% % % % % % % % 
% % % % % % % % 
% % % % % % % % n=length(t1);
% % % % % % % % tt=1:n;
% % % % % % % % nyquist=1/2;
% % % % % % % % freq=(1:n/2)/(n/2)*nyquist;
% % % % % % % % period=1./freq;
% % % % % % % % 
% % % % % % % % for i=1:size(v1,1)
% % % % % % % %     for j=1:size(v1,2)
% % % % % % % %         tg_v1=squeeze(v1(i,j,:));
% % % % % % % %         tg_v2=squeeze(v2(i,j,:));
% % % % % % % %         
% % % % % % % %         fu1=fft(tg_v1);
% % % % % % % %         fu2=fft(tg_v2);
% % % % % % % %         fu1(1) = [];
% % % % % % % %         fu2(1) = [];
% % % % % % % %         power1=abs(fu1(1:floor(n/2))).^2;
% % % % % % % %         power2=abs(fu2(1:floor(n/2))).^2;
% % % % % % % %         
% % % % % % % %         ind1=find(power1==max(power1));
% % % % % % % %         ind2=find(power2==max(power2));
% % % % % % % %         if (period(ind1)==365)
% % % % % % % %             power1(ind1)=0;
% % % % % % % %             ind1=find(power1==max(power1));
% % % % % % % %         end
% % % % % % % %         if (period(ind2)==365)
% % % % % % % %             power2(ind2)=0;
% % % % % % % %             ind2=find(power2==max(power2));
% % % % % % % %         end
% % % % % % % %         if (length(ind1)==0)
% % % % % % % %             p1(i,j)=nan;
% % % % % % % %         else
% % % % % % % %             p1(i,j)=period(ind1);
% % % % % % % %         end
% % % % % % % %         if (length(ind2)==0)
% % % % % % % %             p2(i,j)=nan;
% % % % % % % %         else
% % % % % % % %             p2(i,j)=period(ind2);
% % % % % % % %         end
% % % % % % % %     end
% % % % % % % % end
% % % % % % % %  
% % % % % % % % % p1(find(p1==365))=nan;
% % % % % % % % % p2(find(p2==365))=nan;
% % % % % % % % 
% % % % % % % % figure;hold on;grid on;box on
% % % % % % % % m_proj('mercator','lon',[98 284],'lat',[-20 65]);
% % % % % % % % m_grid('fontsize',20,'fontname','Times','linewidth',2)
% % % % % % % % m_pcolor(lon,lat,p1)
% % % % % % % % h = colorbar;
% % % % % % % % set(h,'fontsize',15)
% % % % % % % % title(h,'day')
% % % % % % % % colormap(jet);
% % % % % % % % % caxis([0 6])
% % % % % % % % m_gshhs_i('color','k')  
% % % % % % % % set(gcf,'Position',[200 100 800 400])
% % % % % % % % saveas(gcf,['D:\figure\sst_period'],'png')
% % % % % % % % 
% % % % % % % % figure;hold on;grid on;box on
% % % % % % % % m_proj('mercator','lon',[98 284],'lat',[-20 65]);
% % % % % % % % m_grid('fontsize',20,'fontname','Times','linewidth',2)
% % % % % % % % m_pcolor(lon,lat,p2)
% % % % % % % % h = colorbar;
% % % % % % % % set(h,'fontsize',15)
% % % % % % % % title(h,'day')
% % % % % % % % colormap(jet);
% % % % % % % % % caxis([0 0.2])
% % % % % % % % m_gshhs_i('color','k')  
% % % % % % % % set(gcf,'Position',[200 100 800 400])
% % % % % % % % saveas(gcf,['D:\figure\std_period'],'png')
% % % % % % % %         
% % % % % % % %   
% % % % % % % %         
% % % % % % % % if (regional_flag ==1)
% % % % % % % %     for box_flag = 1:4
% % % % % % % %         switch box_flag
% % % % % % % %             case 1 % NWP
% % % % % % % %                 xx = [116 154];
% % % % % % % %                 yy = [16 49];
% % % % % % % %                 titlename = 'RMSE NWP (SSHA)';
% % % % % % % %             case 2 % NEP
% % % % % % % %                 xx = [200 280];
% % % % % % % %                 yy = [16 49]
% % % % % % % %                 titlename = 'RMSE NEP (SSHA)';
% % % % % % % %             case 3 % CP
% % % % % % % %                 xx = [160 200];
% % % % % % % %                 yy = [16 49];
% % % % % % % %                 titlename = 'RMSE CP (SSHA)';
% % % % % % % %             case 4 % EQ
% % % % % % % %                 xx = [120 280];
% % % % % % % %                 yy = [-5 5];
% % % % % % % %                 titlename = 'RMSE EQ (SSHA)';
% % % % % % % %             otherwise
% % % % % % % %                 titlename = 'RMSE (SSHA)';
% % % % % % % %         end
% % % % % % % %         ind = find(lon(:,1)>=xx(1) & lon(:,1)<=xx(2) & lat(1,:) >= yy(1) & lat(1,:) <= yy(2));
% % % % % % % %         p1_r(box_flag) = nanmean(p1(ind));
% % % % % % % %         p2_r(box_flag) = nanmean(p2(ind));
% % % % % % % %     end
% % % % % % % %     
% % % % % % % %     bar_x_tick = [1.20:1.1:4.65];
% % % % % % % %     xtick_name = {'NWP'; 'NEP'; 'CP'; 'EQ'};
% % % % % % % %     figure;hold on;grid on;box on
% % % % % % % %     yyaxis left
% % % % % % % %     bar([1:1.1:4.4],p1_r,0.3)
% % % % % % % %     ylabel('day')
% % % % % % % % %     ylim([0 4])
% % % % % % % %     yyaxis right
% % % % % % % %     bar([1.4:1.1:4.9],p2_r,0.3)
% % % % % % % %     ylabel('day')
% % % % % % % % %     ylim([0 0.08])
% % % % % % % %     set(gca,'fontsize',15,'linewidth',2,'xtick',bar_x_tick,'xticklabel',xtick_name)        
% % % % % % % %     legend('period(SST)','period(SSHA)')
% % % % % % % %     saveas(gcf,['D:\figure\period_regional'],'png')
% % % % % % % % 
% % % % % % % % end
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % % % % function [f_final,s_final,phase]=fftspectrum(u,t)
% % % % % % % % 
% % % % % % % % %FFT 
% % % % % % % % % fft(1:10)을 실행하면, 10개의 수가 나오는데, 1번은 0 frequency, 2~5번은 + frequency, 6~10번은
% % % % % % % % % - frequency를 의미한다. 따라서 주파수가 - , 0 , + 순으로 나오게 하려면,
% % % % % % % % % fftshift(fft(1:10))으로 해주면 된다. 0 frequency는 1:10까지의 총합을 의미..
% % % % % % % % 
% % % % % % % % % load tjrwuv2yr.mat
% % % % % % % % 
% % % % % % % % [ns, nt] = size(u);
% % % % % % % % 
% % % % % % % % df = (nt-1)/(max(t)-min(t)); %unit : cycles per day
% % % % % % % % 
% % % % % % % % if mod(nt,2), f=df*[-(nt-1)/2:(nt-1)/2]/nt;
% % % % % % % % else, f = df*[-nt/2:nt/2-1]/nt;
% % % % % % % % end
% % % % % % % % deltaf=diff(f(1:2));
% % % % % % % % 
% % % % % % % % u=u-nanmean(u);      % 0 frequency가 높게나와 그림보는게 어려울 때가 있어 빼준다. 
% % % % % % % % 
% % % % % % % % ii=isnan(u);u(ii)=0; %addhoc
% % % % % % % % %ii = find(isnan(u) | isnan(v));
% % % % % % % % 
% % % % % % % % d=u;
% % % % % % % % fc=fftshift(fft(d))/nt;
% % % % % % % % phase = angle(fc);
% % % % % % % % s=abs(fc).^2/deltaf;
% % % % % % % % 
% % % % % % % % f_final=f(0<f);
% % % % % % % % s_final=s(0<f);
% % % % % % % % %paseval's theorem
% % % % % % % % 
% % % % % % % % % sum(s)*deltaf;  %spetrum의 총 면적(에너지)는 data의 분산값과 일치한다. 단, s에서 평균값을 빼줘야 같다.
% % % % % % % % % nanvar(u);      %분산은 평균값을 포함하지 않기 때문에
% % % % % % % % % 
% % % % % % % % % 
% % % % % % % % 
% % % % % % % % end