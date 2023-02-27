%parameters
N=512;L=8;len=20000;
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

% F is a real-valued function to map the current coe?cient estimate into 
% a certain value of the proportionate step-size parameter
F = @(x)min(2,400*x);
he_final = zeros(length(x),N);
for round=1:10
he = zeros(length(x),N);
for k=N+L+1:len-1
    for i=0:N-1
        Xap(i+1,:)=x(k-i:-1:k+1-L-i);
    end
    hX = he(k,:)*Xap;
    d = h*Xap;
    evector = d-hX+randn(1,L);
    e(k)=mean(evector);
    lInf=max(F(abs(he(k,:))),delta);
    g=max(rou*lInf,F(abs(he(k,:))));
    gm=mean(g);
    g=g/gm;
    G=diag(g);
    he(k+1,:)=he(k,:)+abs((abs(e(k))-abs(e(k-1)))/(abs(e(k))+abs(e(k-1))))*(mu*G*Xap*inv(Xap'*G*Xap+delta*I)*evector')';
end
he_final=he_final+he;
end
he_final=he_final/10;

%%
% figure(1);
% plot(e(N:len).^2,'b');
% xlabel('迭代次数');
% ylabel('均方误差');
% title('PAPA')
% saveas(gcf,'C:\Work\毕设\codes\pic0409\PAPA_MSE.png')


figure(2)
Misalign_EPAPA = 10*log10(sum((he-h).^2,2)/sum(h.^2));
plot(Misalign_EPAPA)
xlabel('迭代次数');
ylabel('系数失调(dB)');
title('EPAPA')
% saveas(gcf,'C:\Work\毕设\codes\pic0409\PAPA_Misalignment.png')





