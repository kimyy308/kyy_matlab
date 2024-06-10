

%ocn

matroot='/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/HCST_EXP/mat/ocn/photoC_TOT_zint_100m';

load([matroot, filesep, 'hcst_pred_ratio_assm_photoC_TOT_zint_100m_v01_v01_l12m.mat']);
pcolor(data2.photoC_TOT_zint_100m_model_ratio_mean_12(:,1:320)'); shading flat; colorbar;

load([matroot, filesep, 'hcst_pred_ratio_assm_photoC_TOT_zint_100m_v01_v01_l13m.mat']);
figure;
pcolor(data2.photoC_TOT_zint_100m_model_ratio_mean_13(:,1:320)'); shading flat; colorbar;


load([matroot, filesep, 'hcst_pred_ratio_assm_photoC_TOT_zint_100m_v01_v01_l37m.mat']);

figure(1); pcolor(data2.photoC_TOT_zint_100m_model_ratio_mean_37(:,1:320)'); shading flat; colorbar; caxis([1 10]);

figure(2); pcolor(abs(data2.photoC_TOT_zint_100m_model_inc_mean_37(:,1:320))'); shading flat; colorbar; caxis([0 0.04]);
figure(3); pcolor(data2.photoC_TOT_zint_100m_model_spr_37(:,1:320)'); shading flat; colorbar; caxis([0 0.04]);

figure(4); pcolor(abs(data2.photoC_TOT_zint_100m_lens2_inc_mean_37(:,1:320))'); shading flat; colorbar; caxis([0 0.04]);
figure(5); pcolor(data2.photoC_TOT_zint_100m_lens2_spr_37(:,1:320)'); shading flat; colorbar; caxis([0 0.04]);


load([matroot, filesep, 'hcst_pred_ratio_assm_photoC_TOT_zint_100m_v01_v01_l36m.mat']);

figure(1); pcolor(data2.photoC_TOT_zint_100m_model_ratio_mean_36(:,1:320)'); shading flat; colorbar; caxis([1 10]);
figure(2); pcolor(data2.photoC_TOT_zint_100m_lens2_ratio_mean_36(:,1:320)'); shading flat; colorbar; caxis([1 10]);
figure(3); pcolor(data2.photoC_TOT_zint_100m_assm_ratio_mean_36(:,1:320)'); shading flat; colorbar; caxis([1 10]);
figure(3); pcolor(mean(data2.photoC_TOT_zint_100m_assm_ratio_36(:,1:320,1:55),3)'); shading flat; colorbar; caxis([1 10]);



figure(2); pcolor(abs(data2.photoC_TOT_zint_100m_model_inc_mean_36(:,1:320))'); shading flat; colorbar; caxis([0 0.04]);
figure(3); pcolor(data2.photoC_TOT_zint_100m_model_spr_36(:,1:320)'); shading flat; colorbar; caxis([0 0.04]);

figure(4); pcolor(abs(data2.photoC_TOT_zint_100m_lens2_inc_mean_36(:,1:320))'); shading flat; colorbar; caxis([0 0.04]);
figure(5); pcolor(data2.photoC_TOT_zint_100m_lens2_spr_36(:,1:320)'); shading flat; colorbar; caxis([0 0.04]);

figure(4); pcolor(abs(data2.photoC_TOT_zint_100m_assm_inc_mean_36(:,1:320))'); shading flat; colorbar; caxis([0 0.04]);
figure(4); pcolor(abs(data2.photoC_TOT_zint_100m_assm_spr_36(:,1:320))'); shading flat; colorbar; caxis([0 0.04]);



%atm
matroot='/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/HCST_EXP/mat/atm/PRECT';

load([matroot, filesep, 'hcst_pred_ratio_assm_PRECT_v01_v01_l37m.mat']);

figure(1); pcolor(data2.PRECT_model_ratio_mean_37(:,1:192)'); shading flat; colorbar; caxis([1 2]);
figure(2); pcolor(data2.PRECT_lens2_ratio_mean_37(:,1:192)'); shading flat; colorbar; caxis([1 2]);
figure(3); pcolor(data2.PRECT_assm_ratio_mean_37(:,1:192)'); shading flat; colorbar; caxis([1 2]);
figure(3); pcolor(mean(data2.PRECT_assm_ratio_37(:,1:192,1:55),3)'); shading flat; colorbar; caxis([1 2]);
