clc; close all; clear all;
load('LIM_Corr.mat')

% ---- data ----
% C  = Corr_ST_Obs{1,1};
% C  = Corr_ST_White{1,1};
% C  = Corr_ST_Colored{1,1};
% C  = Corr_ST_CW{1,1};
C  = Corr_CS_Obs{1,1};
% C  = Corr_CS_White{1,1};     % size: [x y theta lag]
% C  = Corr_CS_Colored{1,1};   % size: [x y theta lag]
% C  = Corr_CS_CW{1,1};          % 예시: 필요에 따라 ST용 변수로 교체하세요 (size: [x y lag]일 수 있음)



fontsize= 20;

ix = 1;                       % NPP index in "current" (column)
targets = [2 3 4 5];          % SST, SSH, MLD, Fe
names   = {'SST','SSH','MLD','Fe'};

% ---- axes & dimension check ----
dims = size(C);
if numel(dims) == 4
    [nX,nY,nTheta,nLag] = size(C);
    hasTheta = true;
else
    [nX,nY,nLag] = size(C);
    nTheta  = 1;
    hasTheta = false;
end
lag = 0:nLag-1;

% ---- figure & layout ----
figure('Color','w');
tlo = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% ---- common settings ----
caxisRange = [-1 1];

if hasTheta
    % ====== CS model: imagesc(θ × lag) 맵 ======
    nColors = 20;                            % discrete steps
    baseColors = [0 0 1; 1 1 1; 1 0 0];      % blue–white–red
    cmap = interp1(linspace(-1,1,3), baseColors, linspace(-1,1,nColors));

    for k = 1:4
        nexttile;
        % M = squeeze(C(ix, targets(k), :, :));      % [theta x lag] (행: θ, 열: lag)
        M = squeeze(C(targets(k), ix, :, :));        % 데이터 구성에 따라 위/아래 중 택1
        imagesc(lag, 1:nTheta, M);
        set(gca,'YDir','normal');
        colormap(gca, cmap);
        caxis(caxisRange);
        title(sprintf('Corr(\\theta, lag): NPP (t+\\tau) \\leftarrow %s (t)', names{k}));
        xlabel('Lead (months)');
        ylabel('\theta (month)');

        if nTheta == 12
            months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
            set(gca,'YTick',1:12,'YTickLabel',months);
        end

        % NaN 투명 처리
        hImg = get(gca,'Children');
        if ~isempty(hImg) && isgraphics(hImg(1),'image')
            set(hImg(1),'AlphaData',~isnan(M));
        end
        grid on;
        set(gca, 'fontsize', fontsize);
    end

    % ---- shared colorbar ----
    cb = colorbar;
    cb.Layout.Tile = 'east';
    cb.Ticks = -1:0.2:1;
    ylabel(cb,'Correlation');
    set(gcf,'Colormap',cmap);

else
    % ====== ST model: lag에 대한 라인 플롯 ======
    for k = 1:4
        nexttile;
        % ST: C의 크기가 [x y lag] 이므로 아래처럼 벡터로 꺼냄
        % (데이터 구성에 따라 (ix, targets(k), :) 또는 (targets(k), ix, :)를 사용)
        yvec = squeeze(C(targets(k), ix, :));  % 필요시 위/아래 전환
        % yvec = squeeze(C(ix, targets(k), :));

        plot(lag, yvec, 'LineWidth', 1.8);
        hold on; yline(0,'k-'); hold off;
        grid on;
        xlim([lag(1) lag(end)]);
        ylim(caxisRange);
        xlabel('Lead (months)');
        ylabel('Correlation');
        title(sprintf('Corr(lag): NPP (t+\\tau) \\leftarrow %s (t)', names{k}));
        set(gca, 'fontsize', fontsize);
    end
    

end

set(gcf, 'Position', [900 500 1200 1000]);

