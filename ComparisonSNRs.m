%parameters
N=512;L=2;len=20000; div=32; M=N/div;
x=randn(1,len);
varepsilon=0.01;
delta=0.01;
rou=0.01;
mu=1;
h=load('randomGenerated_twoCluster_sparsEchopath.mat');
h=h.w0;
h=h';

% X = zeros(1,N);
v = randn(1,length(x));
var_x = zeros(1,length(x));
e = zeros(1,length(x));
e1 = zeros(1,length(x));
I=eye(L);
Xap=zeros(N,L);
h_bs=zeros(div,N/div);
IM=eye(M);
% F is a real-valued function to map the current coe?cient estimate into 
% a certain value of the proportionate step-size parameter
F = @(x)min(2,400*x);

for snr=-2:1
he = zeros(length(x),N);
he1 = zeros(length(x),N);
for k=N+L+1:len
    for i=0:N-1
        Xap(i+1,:)=x(k-i:-1:k+1-L-i);
    end
    
    hX = he(k,:)*Xap;
    hX1 = he1(k,:)*Xap;
    
    h_bs=vec_div(he(k,:),div);
    h_n21=norm21(h_bs);%the norm 21s of the block || row vector
    lInf=max(F(abs(h_n21)),delta);
    g=max(rou*lInf,F(abs(h_n21)));
    gm=mean(g);
    g=g/gm;
    
    h_bs=vec_div(he1(k,:),div);
    h_n21=norm21(h_bs);%the norm 21s of the block || row vector
    lInf=max(F(abs(h_n21)),delta);
    g1=max(rou*lInf,F(abs(h_n21)));
    gm=mean(1);
    g1=g1/gm;
            
    h_G=g(1)*ones(1,M);
    
    for j=2:div
        h_G=[h_G,g(j)*ones(1,M)];%build the block proportionate matrix
    end
    G=diag(h_G);
    
    h_G=g1(1)*ones(1,M);
    
    for j=2:div
        h_G=[h_G,g1(j)*ones(1,M)];%build the block proportionate matrix
    end
    G1=diag(h_G);
    
    d = h*Xap;
    n = 10^snr*randn(1,L);
    evector = d-hX+n;
    e(k)=mean(evector);
    evector1 = d-hX1+n;
    e1(k)=mean(evector1);
    
    he(k+1,:)=he(k,:)+(mu*G*Xap*inv(Xap'*G*Xap+delta*I)*evector')';
    he1(k+1,:)=he1(k,:)+abs((abs(e1(k))-abs(e1(k-1)))/(abs(e1(k))+abs(e1(k-1))))*(mu*G1*Xap*inv(Xap'*G1*Xap+delta*I)*evector1')';
end

Misalign_BSPAPA = 10*log10(sum((he-h).^2,2)/sum(h.^2));
Misalign_BSPAPA1 = 10*log10(sum((he1-h).^2,2)/sum(h.^2));
figure(2)
c=[rand rand rand];
plot(Misalign_BSPAPA,'color',c,'linewidth',1.5);
hold on;
plot(Misalign_BSPAPA1,'color',c,'linewidth',1.5,'LineStyle','--');
hold on;

end
%%
legend('SNR=20dB','SNR=20dB (EDTD)','SNR=10dB','SNR=10dB (EDTD)','SNR=0dB','SNR=0dB (EDTD)','SNR=-10dB','SNR=-10dB (EDTD)')
xlabel('迭代次数');
ylabel('系数失调(dB)');
title('不同信噪比对比')
% saveas(gcf,'C:\Work\毕设\codes\pic0409\PAPA_Misalignment.png')





