F = @(x)min(2,400*x);
%%
N=512;P=8;L=5000;
x=10*randn(1,L);
varepsilon=0.01;
delta=0.01;
rou=0.01;
g=zeros(1,N);
w=load('randomGenerated_twoCluster_sparsEchopath.mat');
w=w.w0;
we = zeros(N,L);
X = zeros(N,P);
%v = randn(1,L);
e = zeros(1,P);
I=eye(P);
MSE=zeros(1,L);
sigma_d=zeros(1,L);
sigma_ye=zeros(1,L);
sigma_v=zeros(1,L);
sigma_e=zeros(P,L);

alphaRECORD = zeros(1,L);
%% P是投影阶数 N是系统长度
K1=0.01;
K2=0.01;
lambda1=1-1/(K1*N);
lambda2=1-1/(K2*N);
amin=0.1;
amax=2;
for t=N+P+1:L
    for i=0:N-1
        X(i+1,:)=x(t-i:-1:t+1-P-i);
    end
    wl=F(abs(we(:,t)));
    lInf=max(wl,delta);
    g=max(rou*lInf,wl);
    gm=mean(g);
    g=g/gm;
    G=diag(g);
    ye=X'*we(:,t);
    d=X'*w;
    e=d-ye+randn(P,1);%iws(2,1.5,1,P)';
    MSE(t)=mean(e).^2;
    sigma_d(t)=lambda2*sigma_d(t-1)+(1-lambda2)*mean(d)^2;
    sigma_ye(t)=lambda2*sigma_ye(t-1)+(1-lambda2)*mean(ye)^2;
    sigma_v(t)=max(0,sigma_d(t)-sigma_ye(t));
    sigma_e(:,t)=(1-lambda1)*sigma_e(:,t-1)+lambda1*e.^2;
    alpha=1-sqrt(sigma_v(t:-1:t+1-P)'./sigma_e(end:-1:1,t));
    alpha=max(alpha,amin);
    alpha=min(alpha,amax);
    alphaRECORD(t)=mean(alpha);
    Alpha=diag(alpha);
    we(:,t+1)=we(:,t)+G*X*inv(X'*G*X+varepsilon*I)*Alpha*e;
end
%%
figure(1)
subplot(211)
plot(w);title('稀疏信道的真实值');
subplot(212)
plot(we(:,L));title('稀疏信道的估计值');
%%
figure(2)
plot(MSE(N:L))
xlabel('迭代次数');
ylabel('均方误差');
title('VSS-SPAPA')
% saveas(gcf,'C:\Work\毕设\codes\pic0409\VSS_SPAPA_MSE.png')


figure(3)
Misalign_VSS_SPAPA = 10*log10(sum((we-w).^2,1)/sum(w.^2));
plot(Misalign_VSS_SPAPA)
xlabel('迭代次数');
ylabel('系数失调(dB)');
title('VSS-SPAPA')
% saveas(gcf,'C:\Work\毕设\codes\pic0409\VSS_SPAPA_Misalignment.png')

figure(4)
stem(1:L,alphaRECORD);
