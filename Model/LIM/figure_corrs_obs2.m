clc; close all; clear all;
load('Obs_LIM_Corr.mat')

% ===== 사용할 상관 맵 선택 =====
% C  = Corr_ST_Obs{51,1};
% C  = Corr_ST_White{51,1};
% C  = Corr_ST_Colored{51,1};
% C  = Corr_ST_CW{51,1};
% C  = Corr_CS_Obs{51,1};      % size: [x y theta lag]  (x=target, y=predictor)
% C  = Corr_CS_White{51,1};
% C  = Corr_CS_Colored{51,1};
C  = Corr_CS_CW{51,1};          % 예시

% ===== 설정 =====
fontsize   = 20;
ix         = 1;                      % NPP index (데이터 정의에 맞게)
targets    = [2 3 4];                % 예: SST, SSH, MLD (Fe 쓰면 바꿔)
names      = {'SST','SSH','MLD'};
caxisRange = [-1 1];
alphaWeak  = 0.25;                   % 비유의 칸 투명도
alphaSig   = 1.00;                   % 유의 칸 투명도
alphaFDR   = 0.1;                   % FDR 유의수준
monthsLbl  = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

% ===== DOF (연도 수) : 1960–2023 → 64년 → 월별 동일 가정 =====
Nyears = 64 * ones(12,1);   % 월별 결측 있으면 여기서 각 월 값으로 수정

% ===== 차원 파악 =====
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

% ===== figure =====
figure('Color','w');
tlo = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% ===== colormap (blue-white-red) =====
nColors    = 20;
baseColors = [0 0 1; 1 1 1; 1 0 0];
cmap       = interp1(linspace(-1,1,3), baseColors, linspace(-1,1,nColors));

if hasTheta
    % ====== CS model: imagesc(θ × lag) ======
    for k = 1:3
        nexttile;
    
        % M 구성 (θ x lag)
        % M = squeeze(C(ix, targets(k), :, :));
        M = squeeze(C(targets(k), ix, :, :));
    
        % 히트맵
        hImg = imagesc(lag, 1:nTheta, M);   % pixel center가 (lag_j, theta_i)
        set(gca,'YDir','reverse');          % Dec 위, Jan 아래
        colormap(gca, cmap); caxis(caxisRange);
        title(sprintf('Corr(\\theta, \\tau): NPP (t+\\tau) \\leftarrow %s (t)', names{k}));
        xlabel('Lead (months)'); ylabel('\theta (month)');
        if nTheta == 12
            set(gca,'YTick',1:12,'YTickLabel',fliplr(monthsLbl));
        end
        grid on; set(gca,'FontSize',fontsize);
    
        % --- 유의성 마스크: DOF 기반 t-검정 + FDR ---
        sigMask = sigmask_from_rmap(M, Nyears, alphaFDR);   % 12 x nLag (true/false)
    
        % --- (새 방식) 유의한 셀에 '점' 찍기 ---
        hold on
        % 1) 부호 구분 없이 모두 '●'로 표시 (간단 버전)
        % [iy, ix] = find(sigMask);           % iy: row(θ), ix: col(τ)
        % scatter(lag(ix), iy, 18, 'k', '.'); % 검은 점, 크기 18
    
        % 2) 옵션: 부호 구분 (양='+', 음='x' 혹은 다른 마커)
        pos = sigMask & (M > 0);
        neg = sigMask & (M < 0);
    
        [iyP, ixP] = find(pos);
        [iyN, ixN] = find(neg);
    
        % 양의 상관: 점('.'), 음의 상관: 'x' 로 표시
        if ~isempty(ixP)
            scatter(lag(ixP), iyP, 20, 'k', '.');   % pos: 굵은 점
        end
        if ~isempty(ixN)
            plot(lag(ixN), iyN, 'kx', 'MarkerSize', 6, 'LineWidth', 1.0); % neg: x 마커
        end
        hold off
    
        % (참고) 알파 마스킹은 사용하지 않음
        % set(hImg, 'AlphaData', ones(size(M)));
    end



    % 공용 컬러바
    cb = colorbar; cb.Layout.Tile = 'east';
    cb.Ticks = -1:0.2:1; ylabel(cb,'Correlation');
    set(gcf,'Colormap',cmap);

else
    % ================================
    % ST: 선그래프 (참고용)
    % ================================
    for k = 1:3
        nexttile;
        yvec = squeeze(C(targets(k), ix, :));
        plot(lag, yvec, 'LineWidth', 1.8);
        hold on; yline(0,'k-'); hold off;
        grid on; xlim([lag(1) lag(end)]); ylim(caxisRange);
        xlabel('Lead (months)'); ylabel('Correlation');
        title(sprintf('Corr(lag): NPP (t+\\tau) \\leftarrow %s (t)', names{k}));
        set(gca,'FontSize',fontsize);
    end
end

set(gcf, 'Position', [900 500 1200 1000]);

% ================================
% 로컬 함수들
% ================================
function sigMask = sigmask_from_rmap(M, Nyears, alphaFDR)
    % M : (12 x nLag) 상관계수 행렬 (월=행, 리드=열)
    % Nyears : (12 x 1) 각 월의 유효 표본수 (여기서는 64로 동일 가정)
    [nTheta, nLag] = size(M);

    % 월별 유효 표본수 행렬 (12 x nLag)
    neff = max(Nyears(:), 3);                 % 최소 3 보장
    neff = repmat(neff(:), 1, nLag);

    % t-통계량 & p값 (양측)
    T = M .* sqrt( (neff - 2) ./ max(1 - M.^2, eps) );
    P = 2 * (1 - tcdf(abs(T), max(neff - 2, 1)));

    % FDR (Benjamini–Hochberg)
    sigMask = fdr_bh_mask(P, alphaFDR);
end

function sigMask = fdr_bh_mask(pvals, alpha)
    p = pvals(:);
    [ps, idx] = sort(p);
    m = numel(p);
    thr = (1:m)'/m * alpha;
    rej = false(m,1);
    rej(idx) = ps <= thr;
    k = find(rej,1,'last');
    rej(:) = false; if ~isempty(k), rej(idx(1:k)) = true; end
    sigMask = reshape(rej, size(pvals));
end
