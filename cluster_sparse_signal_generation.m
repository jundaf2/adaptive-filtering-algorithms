n=512;  
sparsity=0.1;
wpower=32;M1=32/2;

%% random
w0=zeros(n,1);sqrt(wpower/(M1))*randn(n,1);
i=1;
while i<=sparsity*n
posi=fix(rand(1,1)*(n-1))+1;
if(w0(posi)==0)
    w0(posi)=sqrt(wpower/(M1))*randn(1,1);
    i=i+1
end
end
subplot(311)
plot(w0,'k','LineWidth',2);
xlabel('信道长度');
ylabel('幅值');
title('普通（随机）稀疏响应');
%% single cluster
M1=32/2;
M2=16/2;
M3=16;
w=sqrt(wpower/(M1))*zeros(M3,1);
w1=sqrt(wpower/(M2))*randn(M3,1);
A=zeros(256/2,1);
C=zeros(480/2,1);
B=zeros(224/2,1);
w0=[A;C;w;w1;B];
subplot(312)
plot(w0,'k','LineWidth',2);
xlabel('信道长度');
ylabel('幅值');
title('单簇块稀疏响应');
%% two cluster
M1=32/2;
M2=16/2;
M3=16;
w=sqrt(wpower/(M1))*randn(M3,1);
w1=sqrt(wpower/(M2))*randn(M3,1);
A=zeros(256/2,1);
C=zeros(480/2,1);
B=zeros(224/2,1);
w0=[A;w;C;w1;B];
subplot(313)
plot(w0,'k','LineWidth',2);
xlabel('信道长度');
ylabel('幅值');
title('双簇或多簇块稀疏响应');