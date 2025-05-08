function [] = disp_mmspm(mmspm)

disp('=======================')
disp('MMSPM analysis results:')
disp('=======================')
disp(['Low-range slope: ', num2str(-1*mmspm.beta_lo,'%.4f')])
disp(['High-range slope: ', num2str(-1*mmspm.beta_hi,'%.4f')])
disp(['Breakpoint detected at :', num2str(mmspm.breakpoint,'%.2f'), ' Hz'])
disp(['Breakpoint significance: p=', num2str(mmspm.p,'%.4f')])
disp(['Unimodal slope: ', num2str(-1*mmspm.beta_uni)])

end