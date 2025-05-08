function [] = plot_mmspm(mmspm)

p = mmspm.plot_struct;
h = mmspm.p < 0.05;


figure('color','w','units','normalized','outerposition',[0.1,0.15,0.8,0.8]);
subplot(231)
% plot(p.freq,log(p.Y_raw))
hold on
p11 = plot(p.freq,p.Y,'LineWidth',1);
p12 = plot(p.freq,p.E_lo,'LineWidth',1.5);
p13 = plot(p.freq,p.Y_logfit,'LineStyle',':','LineWidth',1.5);
p14 = plot([p.bp_init p.bp_init],[p.Y_logfit(p.dind) p.E_lo(p.dind)],'go-','LineWidth',1.5);
set(gca,'xscale','log','FontSize',12,'LineWidth',1,'box','on')
legend([p11,p12,p13,p14],{'Data','Lower envelope','Dummy fit','Max distance'},'FontSize',12,'location','southwest')
xlabel('Log-Frequency','FontSize',16)
ylabel('Log-Power','FontSize',16)
title('Raw Power Spectrum','FontSize',20)

subplot(232)
p21 = plot(p.freq,p.Yd0,'LineWidth',1);
hold on
p22 = plot(p.freq(p.ind0),p.Yd0(p.ind0),'r*');
set(gca,'FontSize',12,'LineWidth',1,'box','on')
legend([p21,p22],{'Detrended spectrum','Lower quantile'},'FontSize',12,'location','northeast')
xlabel('Frequency','FontSize',16)
ylabel('Log-Power','FontSize',16)
title('Detrended Power Spectrum','FontSize',20)

subplot(233)
plot(p.freq,p.Y,'LineWidth',1)
hold on
p31 = plot(p.f_init,p.Y_init,'r*');
p32 = plot(p.freq,p.A_init,'LineWidth',2);
p33 = plot(mean([p.b0(5) p.b0(5)]),mean([log(p.b0(2)*p.b0(5)^p.b0(1)) log(p.b0(4)*p.b0(5)^p.b0(3))]),'o','color',[0.4 0.4 0.4],'MarkerSize',16,'LineWidth',2);
plot(mean([p.b0(5) p.b0(5)]),mean([log(p.b0(2)*p.b0(5)^p.b0(1)) log(p.b0(4)*p.b0(5)^p.b0(3))]),'.','color',[0.4 0.4 0.4],'MarkerSize',16,'LineWidth',2)
set(gca,'xscale','log','FontSize',12,'LineWidth',1,'box','on')
legend([p31,p32,p33],{'Lower quantile','Initial bipolar fit', 'Breakpoint estimate'},'FontSize',12,'location','southwest')
xlabel('Log-Frequency','FontSize',16)
ylabel('Log-Power','FontSize',16)
title('Initial Bimodal Fit','FontSize',20)

subplot(234)
p41 = plot(p.freq,p.Yd,'LineWidth',1);
hold on
p42 = plot(p.freq,p.gf_curve,'LineWidth',2);
set(gca,'FontSize',12,'LineWidth',1,'box','on')
legend([p41,p42],{'Data','Multi-Gaussian fit'},'FontSize',12,'location','northeast')
xlabel('Frequency','FontSize',16)
ylabel('Log-Power','FontSize',16)
title('Oscillatory Peaks','FontSize',20)

subplot(235)
p51 = plot(p.freq,p.Y,'LineWidth',1);
hold on
% p52 = plot(freq_np,Y_np,'r*');
p52 = plot(p.freq,p.Yf,'LineWidth',1.5);
set(gca,'xscale','log','FontSize',12,'LineWidth',1,'box','on')
legend([p51,p52],{'Data','Peaks removed'},'FontSize',12,'location','southwest')
xlabel('Log-Frequency','FontSize',16)
ylabel('Log-Power','FontSize',16)
title('Aperiodic Spectrum','FontSize',20)

subplot(236)
p61 = plot(p.freq,p.Y,'LineWidth',1);
hold on
p62 = plot(p.freq,p.Y_mmspm,'LineWidth',1.5);
p63 = plot(p.freq,p.A,'LineWidth',2);
if h
    p64 = plot(mean([p.b(5) p.b(5)]),mean([log(p.b(2)*p.b(5)^p.b(1)) log(p.b(4)*p.b(5)^p.b(3))]),'go','MarkerSize',16,'LineWidth',2);
    plot(mean([p.b(5) p.b(5)]),mean([log(p.b(2)*p.b(5)^p.b(1)) log(p.b(4)*p.b(5)^p.b(3))]),'g.','MarkerSize',16,'LineWidth',2)
else
    p64 = plot(mean([p.b(5) p.b(5)]),mean([log(p.b(2)*p.b(5)^p.b(1)) log(p.b(4)*p.b(5)^p.b(3))]),'o','color',[0.4 0.4 0.4],'MarkerSize',16,'LineWidth',2);
    plot(mean([p.b(5) p.b(5)]),mean([log(p.b(2)*p.b(5)^p.b(1)) log(p.b(4)*p.b(5)^p.b(3))]),'.','color',[0.4 0.4 0.4],'MarkerSize',16,'LineWidth',2)
end
set(gca,'xscale','log','FontSize',12,'LineWidth',1,'box','on')
legend([p61,p62,p63,p64],{'Data','MM-SPM spectrum','Bipolar fit', 'Breakpoint'},'FontSize',12,'location','southwest')
xlabel('Log-Frequency','FontSize',16)
ylabel('Log-Power','FontSize',16)
title('MM-SPM & Final Bimodal Fit','FontSize',20)


end