%'NLMS','PNLMS','IPNLMS'
%parameters
N=512;L=30000;
x=10*randn(1,L);
varepsilon=0.01;
sigma=0.01;
delta=0.01;
rou=5/N;
theta=0;
mu=0.5;
h=load('randomGenerated_twoCluster_sparsEchopath.mat');
h=h.w0;
h=h';
h=zeros(1,length(h));
spst=10;
h(100:100-1+spst)=randn(1,spst);


X = zeros(1,N);
v = randn(1,length(x));
e = zeros(1,length(x));
var_x = zeros(1,length(x));
g=zeros(1,N);
y=zeros(1,L);
ye=zeros(1,L);

he_NLMS = zeros(length(x),N);
he_PNLMS = zeros(length(x),N);
he_IPNLMS = zeros(length(x),N);

he_NLMS_final = zeros(length(x),N);
he_PNLMS_final = zeros(length(x),N);
he_IPNLMS_final = zeros(length(x),N);
for round=1:1
for k=N+1:L-1
    noise=randn;
    X = x(k:-1:k-N+1);
    
    y(k)= h*X';
    
    var_x(k) = 1/N*sum(X.^2);
    
    %PNLMS
    ye(k)= he_PNLMS(k,:)*X';
    e(k)=y(k)-ye(k)+noise;
    lInf=max(max(abs(he_PNLMS(k,:))),delta);
    g=max(rou*lInf,abs(he_PNLMS(k,:)));
    gm=mean(g);
    he_PNLMS(k+1,:)=he_PNLMS(k,:)+mu/N/gm*e(k)/var_x(k)*g.*X;
    
    %IPNLMS
    ye(k)= he_IPNLMS(k,:)*X';
    e(k)=y(k)-ye(k)+noise;
    lInf=max(max(abs(he_IPNLMS(k,:))),delta);
     g = (1-theta)/(2*L)+(1+theta)*abs(he_IPNLMS(k,:))/(2*norm(he_IPNLMS(k,:),1)+varepsilon);
    gm=mean(g);
    he_IPNLMS(k+1,:)=he_IPNLMS(k,:)+mu/N/gm*e(k)/var_x(k)*g.*X;
    
    %NLMS
    ye(k)= he_NLMS(k,:)*X';
    e(k)=y(k)-ye(k)+noise;
    he_NLMS(k+1,:)=he_NLMS(k,:)+mu/N*e(k)/var_x(k)*X;
end
he_NLMS_final=he_NLMS_final+he_NLMS;
he_PNLMS_final=he_PNLMS_final+he_PNLMS;
he_IPNLMS_final=he_IPNLMS_final+he_IPNLMS;
end
he_NLMS_final=he_NLMS_final/round;
he_PNLMS_final=he_PNLMS_final/round;
he_IPNLMS_final=he_IPNLMS_final/round;

%%
% figure(1);
% plot(e(N:L).^2,'b');
% xlabel('迭代次数');
% ylabel('均方误差');
% title('NLMS')
% saveas(gcf,'C:\Work\毕设\codes\pic0409\NLMS_MSE.png')


figure(2)
Misalign_NLMS = 10*log10(sum((he_NLMS_final-h).^2,2)/sum(h.^2));
Misalign_PNLMS = 10*log10(sum((he_PNLMS_final-h).^2,2)/sum(h.^2));
Misalign_IPNLMS = 10*log10(sum((he_IPNLMS_final-h).^2,2)/sum(h.^2));
plot(Misalign_NLMS,'k','linewidth',1.5)
hold on
plot(Misalign_PNLMS,'b','linewidth',1.5)
hold on
plot(Misalign_IPNLMS,'r','linewidth',1)
hold on
legend('NLMS','PNLMS','IPNLMS')
title('NMSE')
xlabel('迭代次数');
ylabel('归一化系数失调(dB)');
% saveas(gcf,'C:\Work\毕设\codes\pic0409\NLMS_Misalignment.png')


