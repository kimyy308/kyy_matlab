clc; close all; clear all;
load('/Volumes/kyy_raid/kimyy/Research/Postdoc/03_IBS/2025_LIM_BGC_prediction/LIM_Corr.mat')

% ===== 사용할 상관 맵 선택 =====
% C  = Corr_ST_Obs{1,1};
% C  = Corr_ST_White{1,1};
% C  = Corr_ST_Colored{1,1};
% C  = Corr_ST_CW{1,1};
% C  = Corr_CS_Obs{1,1};         % size: [x y theta lag]  (x=target, y=predictor)
% C  = Corr_CS_White{1,1};
% C  = Corr_CS_Colored{1,1};
C  = Corr_CS_CW{1,1};

% ===== 설정 =====
fontsize   = 20;
ix         = 1;                        % NPP index (데이터 정의에 맞게)
targets    = [2 3 4 5];                % SST, SSH, MLD, Fe
names      = {'SST','SSH','MLD','Fe'};
titles     = {'(a)', '(b)', '(c)', '(d)'};
caxisRange = [-1 1];
alphaFDR   = 0.05;                     % FDR 유의수준 (q)
monthsLbl  = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

% ===== DOF(표본 수) 설정 =====
% 예: 1960–2023 → 64년
Nyears = 64 * ones(12,1);

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

% ===== figure & colormap =====
figure('Color','w');
tlo = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nColors    = 20;
baseColors = [0 0 1; 1 1 1; 1 0 0];
cmap       = interp1(linspace(-1,1,3), baseColors, linspace(-1,1,nColors));

if hasTheta
    % ====== CS model: imagesc(θ × lag) ======
    for k = 1:4
        nexttile;

        M = squeeze(C(ix, targets(k), :, :));      % [theta x lag]
%         M = squeeze(C(targets(k), ix, :, :));        % [theta x lag]

        % 히트맵
        imagesc(lag, 1:nTheta, M);
        set(gca,'YDir','normal');                    % Jan 위, Dec 아래
        colormap(gca, cmap); caxis(caxisRange);
        title(sprintf('%s K(\\theta, \\tau): NPP (t+\\tau) \\leftarrow %s (t)', titles{k}, names{k}));

        xlabel('Lead (months)'); ylabel('\theta (month)');
        if nTheta == 12
            set(gca,'YTick',1:12,'YTickLabel',monthsLbl);
        end
        grid on; set(gca,'FontSize',fontsize);

%         % --- 유의성: DOF 기반 t-검정 + FDR (원데이터 불필요) ---
%         sigMask = sigmask_from_rmap(M, Nyears, alphaFDR);   % 12 x nLag (true/false)

        % fixed threshold
        N = 64; nu = N-2;
        tcrit = tinv(1-0.05/2, nu);                    % 양측 0.05 예시
        rcrit = sqrt( tcrit^2 / (tcrit^2 + nu) );      % ≈ 0.246 (N=64)
        sigMask = abs(M) >= rcrit;

        % --- 유의 셀에 동일한 '점'만 찍기 (부호 구분 없음) ---
        [iy, ixDot] = find(sigMask);                 % iy=row(θ), ixDot=col(τ)
        hold on
        if ~isempty(ixDot)
            scatter(lag(ixDot), iy, 40, 'k', '.');   % 검정 점(●), 크기 20
            % 필요시 색/크기: scatter(lag(ixDot), iy, 18, 'w', '.');  % 흰 점
        end
        hold off
    end

    % 공용 컬러바
    cb = colorbar; cb.Layout.Tile = 'east';
    cb.Ticks = -1:0.2:1; ylabel(cb,'Correlation');
    set(gcf,'Colormap',cmap);

else
    % ====== ST model: lag 선그래프 (참고) ======
    for k = 1:4
        nexttile;
        yvec = squeeze(C(targets(k), ix, :));
        plot(lag, yvec, 'LineWidth', 1.8);
        hold on; yline(0,'k-'); hold off;
        grid on; xlim([lag(1) lag(end)]); ylim(caxisRange);
        xlabel('Lead (months)'); ylabel('Correlation');
        title(sprintf('K(lag): NPP (t+\\tau) \\leftarrow %s (t)', names{k}));
        set(gca,'FontSize',fontsize);
    end
end

set(gcf, 'Position', [900 500 1200 1000]);

% ================================
% 로컬 함수 (유의성)
% ================================
function sigMask = sigmask_from_rmap(M, Nyears, alphaFDR)
    % M : (12 x nLag) 상관계수 행렬 (월=행, 리드=열)
    % Nyears : (12 x 1) 각 월의 유효 표본수
    [nTheta, nLag] = size(M);

    neff = max(Nyears(:), 3);
    neff = repmat(neff(:), 1, nLag);

    T = M .* sqrt( (neff - 2) ./ max(1 - M.^2, eps) );
    P = 2 * (1 - tcdf(abs(T), max(neff - 2, 1)));

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
