pd1 = makedist('Stable','alpha',2,'beta',0,'gam',1,'delta',0);
pd2 = makedist('Stable','alpha',1,'beta',0,'gam',1,'delta',0);
pd3 = makedist('Stable','alpha',0.5,'beta',0,'gam',1,'delta',0);

x = -5:.1:5;
pdf1 = pdf(pd1,x);
pdf2 = pdf(pd2,x);
pdf3 = pdf(pd3,x);

figure
plot(x,pdf1,'b-','linewidth',2);
hold on
plot(x,pdf2,'r-.','linewidth',2);
plot(x,pdf3,'k--','linewidth',2);
title('不同\alpha参数下的稳定分布概率密度函数')
legend('\alpha = 2','\alpha = 1','\alpha = 0.5')
hold off