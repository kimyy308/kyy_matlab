clear all; close all; clc;


x=[1:100].*0.5;
plot(x, cos(x), 'k', 'linewidth', 2);
hold on;
plot(x, [cos(x(1)), cos(x(2:100)).*0.6], 'b', 'linewidth', 2);
% yline(cos(0.5), 'r', 'linewidth', 2);
plot(x, repmat(cos(x(1)), 1, 100), 'r', 'linewidth', 2)
hold off
xlabel('Month')
ylabel('Value')
legend({'OBS', 'Model A', 'Model B'}, 'Location', 'SouthOutside', 'Orientation', 'horizontal' )
set(gca, 'Fontsize', 15)

AR1coef=corrcoef(cos(x(1:99)), cos(x(2:100)));
AR1(1)=cos(x(1));
for i=2:100
    AR1(i)=AR1(i-1).*AR1coef(1,2);
end

x=[1:100].*0.5;
plot(x, cos(x), 'k', 'linewidth', 2);
hold on;
plot(x, [cos(x(1)), cos(x(2:100)).*0.6], 'b', 'linewidth', 2);
% yline(cos(0.5), 'r', 'linewidth', 2);
plot(x, AR1, 'r', 'linewidth', 2)
hold off
xlabel('Month')
ylabel('Value')
legend({'OBS', 'Model A', 'Model B'}, 'Location', 'SouthOutside', 'Orientation', 'horizontal' )
set(gca, 'Fontsize', 15)


%% initial error case (mean biased)
x=[1:100].*0.5;
plot(x, cos(x), 'k', 'linewidth', 2);
hold on;
plot(x, AR1, 'r', 'linewidth', 2)
plot(x, [cos(x(1)), cos(x(2:100)).*0.6]-1, 'b', 'linewidth', 2);
% yline(cos(0.5), 'r', 'linewidth', 2);
% plot(x, repmat(cos(x(1)), 1, 100), 'r', 'linewidth', 2)
hold off
xlabel('Month')
ylabel('Value')
legend({'OBS', 'Model B', 'Model C'}, 'Location', 'SouthOutside', 'Orientation', 'horizontal' )
set(gca, 'Fontsize', 15)


%% initial error case (error in phasing)
x=[1:100].*0.5;
noise=(rand(1,length(AR1))*2-1)*0.05;
AR1_noise=AR1+noise;
plot(x, cos(x), 'k', 'linewidth', 2);
hold on;
plot(x, AR1_noise, 'r', 'linewidth', 2)
plot(x, cos(x+2), 'b', 'linewidth', 2);
% yline(cos(0.5), 'r', 'linewidth', 2);
% plot(x, repmat(cos(x(1)), 1, 100), 'r', 'linewidth', 2)
hold off
xlabel('Month')
ylabel('Value')
% legend({'OBS', 'AR1', 'Hindcast'}, 'Location', 'SouthOutside', 'Orientation', 'horizontal' )
legend({'OBS', 'Model B', 'Model C'}, 'Location', 'SouthOutside', 'Orientation', 'horizontal' )

set(gca, 'Fontsize', 15)

sqrt(sum((cos(x)-AR1_noise).^2))
sqrt(sum((cos(x)-cos(x+2)).^2))


%% uncertainty growth in initialized ensemble forecast
x=[1:100].*0.5;

for j=1:100 % ensemble member
    y(1,j)=0.5;
    z(1,j)=-1+2*pi*j/100.0;
    for i=2:100 % time
        z(i,j)=z(i-1,j)+0.5+(2*rand(1)-1)*0.7;
        y(i,j)=y(i-1,j)+0.5+(2*rand(1)-1)*0.7;
    end
end
z_em=mean(cos(z),2);
z_es=1*std(cos(z),1,2);
z_mx=max(cos(z),[],2);
z_mn=min(cos(z),[],2);
z_low=z_em-z_es/2; z_upper=z_em+z_es/2; %std based
% z_low=z_mn; z_upper=z_mx/2; %minmax based

y_em=mean(cos(y),2);
y_es=1*std(cos(y),1,2);
y_mx=max(cos(y),[],2);
y_mn=min(cos(y),[],2);
y_low=y_em-y_es/2; y_upper=y_em+y_es/2; %std based
% y_low=y_mn; y_upper=y_mx; %minmax based


%% colorset
val_transparent = 0.15;
cmap_HCST = [0,0,1]; % blue
cmap_HCST_b = rgb2hsv(cmap_HCST);
cmap_HCST_b(:,2) =  val_transparent;
cmap_HCST_b = hsv2rgb(cmap_HCST_b);

cmap_ASSM = [1, 0, 0]; % Red [1, 0, 0], Yellow [1, 1, 0]
cmap_ASSM_b = rgb2hsv(cmap_ASSM);
cmap_ASSM_b(:,2) =  val_transparent;
cmap_ASSM_b = hsv2rgb(cmap_ASSM_b);

cmap_LENS2 = [0.5, 1.0, 0.5]; % Red [1, 0, 0], Yellow [1, 1, 0], Green [0.5, 1.0, 0.5]
cmap_LENS2_b = rgb2hsv(cmap_LENS2);
cmap_LENS2_b(:,2) =  val_transparent;
cmap_LENS2_b = hsv2rgb(cmap_LENS2_b);

%% initialized
% spread
fig_ts.range=fill([x, flip(x)], ...
                [y_low(1:end); flip(y_upper(1:end))], cmap_HCST_b);
fig_ts.range.EdgeColor = 'none';
set(get(get(fig_ts.range,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

hold on
% observation
plot(x, cos(x), 'k', 'linewidth', 2);

% hcst ensemble mean
plot(x, y_em, 'color', cmap_HCST, 'linewidth', 2);
xlabel('Month')
ylabel('Value')
legend({'OBS', 'Initialized'}, 'Location', 'SouthOutside', 'Orientation', 'horizontal' )


% dots of ens members
for j=1:100
    fig_sc=scatter(x(1), cos(y(1,j)), 'b*');
    set(get(get(fig_sc,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
for j=1:100
    fig_sc=scatter(x(80), cos(y(80,j)), 'b*');
    set(get(get(fig_sc,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end

hold off
set(gca, 'Fontsize', 15)




%% uninitialized
% spread
fig_ts.range=fill([x, flip(x)], ...
                [z_low(1:end); flip(z_upper(1:end))], cmap_ASSM_b);
fig_ts.range.EdgeColor = 'none';
set(get(get(fig_ts.range,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

hold on
% observation
plot(x, cos(x), 'k', 'linewidth', 2);

% hcst ensemble mean
plot(x, z_em, 'color', cmap_ASSM, 'linewidth', 2);

% dots of ens members
for j=1:100
    fig_sc=scatter(x(1), cos(z(1,j)), 'r*');
    set(get(get(fig_sc,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
for j=1:100
    fig_sc=scatter(x(80), cos(z(80,j)), 'r*');
    set(get(get(fig_sc,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end

hold off
xlabel('Month')
ylabel('Value')
legend({'OBS', 'Uninitialized'}, 'Location', 'SouthOutside', 'Orientation', 'horizontal' )
set(gca, 'Fontsize', 15)




plot(x, cos(x), 'k', 'linewidth', 2);
hold on;
plot(x, cos(z), 'r', 'linewidth', 2)
plot(x, cos(y), 'b', 'linewidth', 2);
% yline(cos(0.5), 'r', 'linewidth', 2);
% plot(x, repmat(cos(x(1)), 1, 100), 'r', 'linewidth', 2)
hold off
xlabel('Month')
ylabel('Value')
legend({'OBS', 'Uninitilized', 'Initialized'}, 'Location', 'SouthOutside', 'Orientation', 'horizontal' )
set(gca, 'Fontsize', 15)



ti=50;
x=[1:ti].*0.5;
plot(x, cos(x)+x, 'k', 'linewidth', 2);
hold on;
% plot(x, AR1+x, 'r', 'linewidth', 2)
plot(x, [cos(x(1)), cos(x(2:ti)).*0.6]-1+x, 'b', 'linewidth', 2);
% yline(cos(0.5), 'r', 'linewidth', 2);
% plot(x, repmat(cos(x(1)), 1, 100), 'r', 'linewidth', 2)
hold off
xlabel('Year')
ylabel('Value')
legend({'OBS', 'Model D'}, 'Location', 'SouthOutside', 'Orientation', 'horizontal' )
set(gca, 'Fontsize', 15)



%% random noise
y=rand(1,50);
plot(x,[y(1:47), NaN(1,3)],'linewidth', 2)
xlabel('Month')
ylabel('Value')
legend({'Variable A'}, 'Location', 'SouthOutside', 'Orientation', 'horizontal' )
set(gca, 'Fontsize', 15)

%% signal
y=cos(x/3);
plot(x,[y(1:47), NaN(1,3)],'linewidth', 2)
xlabel('Month')
ylabel('Value')
legend({'Variable A'}, 'Location', 'SouthOutside', 'Orientation', 'horizontal' )
set(gca, 'Fontsize', 15)


%% signal + noise
y=cos(x/3) + (rand(1,50)*2-1);
plot(x,[y(1:47), NaN(1,3)],'linewidth', 2)
xlabel('Month')
ylabel('Value')
legend({'Variable A'}, 'Location', 'SouthOutside', 'Orientation', 'horizontal' )
set(gca, 'Fontsize', 15)