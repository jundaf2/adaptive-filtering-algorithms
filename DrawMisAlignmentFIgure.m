
m1=load('Misalign_NLMS');
m2=load('Misalign_PAPA');
m3=load('Misalign_PNLMS');
m4=load('Misalign_VSS_SPAPA');
m5=load('Misalign_VSS_IPAPSA');
m6=load('Misalign_EPAPA');
figure();
% plot(m1.Misalign_NLMS,'r');hold on;
plot(m2.Misalign_PAPA,'g');hold on;
% plot(m3.Misalign_PNLMS,'b');hold on;
% plot(m4.Misalign_VSS_SPAPA,'k');hold on;
% plot(m5.Misalign_VSS_IPAPSA,'y');hold on;
plot(m6.Misalign_EPAPA,'b');hold on;


xlabel('迭代次数');
ylabel('系数失调(dB)');
% title('NLMS vs PNLMS vs PAPA vs VSS-SPAPA vs VSS-IPAPSA')
% legend('NLMS','PAPA','PNLMS','VSS-SPAPA','VSS-IPAPSA')

