%parameters
N=512;L=4;len=20000;
x=10*randn(1,len);
varepsilon=0.01;
delta=0.01;
rou=0.01;
mu=1;
g(N)=0;
h=load('randomGenerated_twoCluster_sparsEchopath.mat');
h=h.w0;
h=h';

% X = zeros(1,N);
v = randn(1,length(x));
var_x = zeros(1,length(x));
e = zeros(1,length(x));
I=eye(L);
Xap=zeros(N,L);

h=zeros(1,length(h));
h(100:109)=randn(1,10);
%%
he_NLMS = zeros(length(x),N);
X = zeros(1,N);
v = randn(1,length(x));
var_x = zeros(1,length(x));
e = zeros(1,length(x));
y=zeros(1,L);
ye=zeros(1,L);
he_NLMS_final = zeros(length(x),N);
for round=1:10
for k=N+1:L-1
    X = x(k:-1:k-N+1);
    ye(k)= he(k,:)*X';
    y(k)= h*X';
    e(k)=y(k)-ye(k)+randn;
    var_x(k) = 1/N*sum(X.^2);
    
    he_NLMS(k+1,:)=he_NLMS(k,:)+mu/N*e(k)/var_x(k)*X;
end
he_NLMS_final=he_NLMS_final+he_NLMS;
end
he_NLMS_final=he_NLMS_final/round;

%%
he_APA_final = zeros(length(x),N);
he_APA = zeros(length(x),N);
for round=1:10
for k=N+L+1:len-1
    for i=0:N-1
        Xap(i+1,:)=x(k-i:-1:k+1-L-i);
    end
    hX = he(k,:)*Xap;
    d = h*Xap;
    evector = d-hX+randn(1,L);
    e(k)=mean(evector);

    he_APA(k+1,:)=he_APA(k,:)+(mu*Xap*inv(Xap'*Xap+delta*I)*evector')';
end
he_APA_final=he_APA_final+he_APA;
end
he_APA_final=he_APA_final/round;

%%
% figure(1);
% plot(e(N:len).^2,'b');
% xlabel('迭代次数');
% ylabel('均方误差');
% title('APA')
% saveas(gcf,'C:\Work\毕设\codes\pic0409\APA_MSE.png')


figure(2)
Misalign_APA = 10*log10(sum((he_final-h).^2,2)/sum(h.^2));
plot(Misalign_APA)
hold on
xlabel('迭代次数');
ylabel('系数失调(dB)');
title('APA')
% saveas(gcf,'C:\Work\毕设\codes\pic0409\APA_Misalignment.png')





