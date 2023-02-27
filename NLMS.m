%parameters
N=512;L=30000;
x=10*randn(1,L);
varepsilon=0.01;
sigma=0.01;
mu=0.5;
h=load('randomGenerated_twoCluster_sparsEchopath.mat');
h=h.w0;
h=h';
h=zeros(1,length(h));
h(100:109)=randn(1,10);
he = zeros(length(x),N);
X = zeros(1,N);
v = randn(1,length(x));
var_x = zeros(1,length(x));
e = zeros(1,length(x));
y=zeros(1,L);
ye=zeros(1,L);
he_final = zeros(length(x),N);
for round=1:10
for k=N+1:L-1
    X = x(k:-1:k-N+1);
    ye(k)= he(k,:)*X';
    y(k)= h*X';
    e(k)=y(k)-ye(k)+randn;
    var_x(k) = 1/N*sum(X.^2);
    
    
    he(k+1,:)=he(k,:)+mu/N*e(k)/var_x(k)*X;
end
he_final=he_final+he;
end
he_final=he_final/round;
%%
% figure(1);
% plot(e(N:L).^2,'b');
% xlabel('迭代次数');
% ylabel('均方误差');
% title('NLMS')
% saveas(gcf,'C:\Work\毕设\codes\pic0409\NLMS_MSE.png')


figure(2)
Misalign_NLMS = 10*log10(sum((he_final-h).^2,2)/sum(h.^2));
plot(Misalign_NLMS,'k')
legend('NLMS')
hold on
xlabel('迭代次数');
ylabel('系数失调(dB)');
title('NLMS')
saveas(gcf,'C:\Work\毕设\codes\pic0409\NLMS_Misalignment.png')


