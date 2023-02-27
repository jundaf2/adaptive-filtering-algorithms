%parameters
N=512;L=8;len=20000;
x=10*randn(1,len);
%%
d=[1 -1.5 0.7 0.1]; c=[1 0.5 0.2];  % 分子分母多项式系数
nd=length(d)-1 ;nc=length(c)-1;  %阶次
xik=zeros(nc,1);  %白噪声初值
ek=zeros(nd,1);
xi=randn(len,1);  %产生均值为0，方差为1的高斯白噪声序列

%%
varepsilon=0.01;
delta=0.01;
rou=0.01;
mu=1;
g(N)=0;
h=load('randomGenerated_twoCluster_sparsEchopath.mat');
h=h.w0;
h=h';

% X = zeros(1,N);
ROUND = 20;
I=eye(L);
Xap=zeros(N,L);

h=zeros(1,length(h));
h(100:109)=randn(1,10);
%% NLMS
he_NLMS = zeros(length(x),N);
X = zeros(1,N);
v = randn(1,length(x));
var_x = zeros(1,length(x));
e = zeros(1,length(x));
y=zeros(1,len);
ye=zeros(1,len);
he_NLMS_final = zeros(length(x),N);
for round=1:ROUND
for k=N+1:len-1
    X = x(k:-1:k-N+1);
    ye(k)= he_NLMS(k,:)*X';
    y(k)= h*X';
    e(k)=y(k)-ye(k)+randn;
    var_x(k) = 1/N*sum(X.^2);
    
    he_NLMS(k+1,:)=he_NLMS(k,:)+mu/N*e(k)/var_x(k)*X;
end
he_NLMS_final=he_NLMS_final+he_NLMS;
end
he_NLMS_final=he_NLMS_final/ROUND;

%% EDTD-NLMS
he_NLMS1 = zeros(length(x),N);
X = zeros(1,N);
v = randn(1,length(x));
var_x = zeros(1,length(x));
e = zeros(1,length(x));
y=zeros(1,len);
ye=zeros(1,len);
he_NLMS_final1 = zeros(length(x),N);
for round=1:ROUND
for k=N+1:len-1
    X = x(k:-1:k-N+1);
    ye(k)= he_NLMS1(k,:)*X';
    y(k)= h*X';
    e(k)=y(k)-ye(k)+randn;
    var_x(k) = 1/N*sum(X.^2);
    
    he_NLMS1(k+1,:)=he_NLMS1(k,:)+abs(e(k))-abs(e(k-1))*mu/N*e(k)/var_x(k)*X;
end
he_NLMS_final1=he_NLMS_final1+he_NLMS1;
end
he_NLMS_final1=he_NLMS_final1/ROUND;

%%


figure(2)

Misalign_NLMS = 10*log10(sum((he_NLMS_final-h).^2,2)/sum(h.^2));

Misalign_APA = 10*log10(sum((he_NLMS_final1-h).^2,2)/sum(h.^2));

plot(Misalign_NLMS,'k','linewidth',1.5)
hold on
plot(Misalign_APA,'b','linewidth',1.5)
hold on

legend('NLMS','EDTD-NLMS')
title('NMSE')
xlabel('迭代次数');
ylabel('归一化系数失调(dB)');
