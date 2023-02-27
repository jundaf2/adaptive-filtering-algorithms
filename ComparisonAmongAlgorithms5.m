% EDTD-BS-PAPSA 的 跟踪性能
N=512;L=2;len=200000; 
x=5*randn(1,len);
% a1=-0.5;
% for k=2:len
%      x(k)=-a1*x(k-1)+x(k);
% end 

varepsilon=0.01;
delta=0.01;
rou=0.01;
mu=1;
h=load('randomGenerated_twoCluster_sparsEchopath.mat');
h=h.w0;
h=h';

figure (1);
subplot(211)


h=zeros(1,length(h));
spst=12;
h(128:128-1+spst)=5*randn(1,spst);
h(256:256-1+spst)=5*randn(1,spst);
plot(h,'k','linewidth',1.5);
xlabel('(a) 双簇信道')

subplot(212)
h1=zeros(1,length(h));
spst=24;
h1(342:342-1+spst)=5*randn(1,spst);
plot(h1,'k','linewidth',1.5);
xlabel('(b) 单簇信道')

%%
ROUND=10;
% X = zeros(1,N);
% v = randn(1,length(x));
% var_x = zeros(1,length(x));
% e = zeros(1,length(x));
I=eye(L);
Xap=zeros(N,L);

% F is a real-valued function to map the current coe?cient estimate into 
% a certain value of the proportionate step-size parameter
F = @(x)min(2,400*x);

a = 1;gama = 1;ms = 1;
vn = 2*iws(a,gama,ms,len);


%% BSPAPA 
div=16; M=N/div;
he_BSPAPA1 = zeros(length(x),N);
h_bs=zeros(div,M);
IM=eye(M);


he_BSPAPA1_final = zeros(length(x),N);
for round=1:ROUND
round
for k=N+L+1:len/2-1
    for i=0:N-1
        Xap(i+1,:)=x(k-i:-1:k+1-L-i);
    end
    
    hX = he_BSPAPA1(k,:)*Xap;
    
    h_bs=vec_div(he_BSPAPA1(k,:),div);
    h_n21=norm21(h_bs);%the norm 21s of the block || row vector
    lInf=max(abs(h_n21),delta);
    g=max(rou*lInf,abs(h_n21));
    gm=mean(g);
    g=g/gm;
            
    h_G=g(1)*ones(1,M);
    
    for j=2:div
        h_G=[h_G,g(j)*ones(1,M)];%build the block proportionate matrix
    end
    G=diag(h_G);
    
    d = h*Xap;
    evector = d-hX+vn(k:-1:k-L+1);
    
    he_BSPAPA1(k+1,:)=he_BSPAPA1(k,:)+(mu*G*Xap*inv(Xap'*G*Xap+delta*I)*sign(evector'))';
end

for k=len/2:len-1
    for i=0:N-1
        Xap(i+1,:)=x(k-i:-1:k+1-L-i);
    end
    
    hX = he_BSPAPA1(k,:)*Xap;
    
    h_bs=vec_div(he_BSPAPA1(k,:),div);
    h_n21=norm21(h_bs);%the norm 21s of the block || row vector
    lInf=max(abs(h_n21),delta);
    g=max(rou*lInf,abs(h_n21));
    gm=mean(g);
    g=g/gm;
            
    h_G=g(1)*ones(1,M);
    
    for j=2:div
        h_G=[h_G,g(j)*ones(1,M)];%build the block proportionate matrix
    end
    G=diag(h_G);
    
    d = h1*Xap;
    evector = d-hX+vn(k:-1:k-L+1);
    
    he_BSPAPA1(k+1,:)=he_BSPAPA1(k,:)+(mu*G*Xap*inv(Xap'*G*Xap+delta*I)*sign(evector'))';
end

he_BSPAPA1_final=he_BSPAPA1_final+he_BSPAPA1;
end
he_BSPAPA1_final=he_BSPAPA1_final/ROUND;


%% EDTD BSPAPA 
div=16; M=N/div;
he_BSPAPSA = zeros(length(x),N);
h_bs=zeros(div,M);
IM=eye(M);


he_BSPAPSA_final = zeros(length(x),N);
for round=1:ROUND
round
for k=N+L+1:len/2-1
    for i=0:N-1
        Xap(i+1,:)=x(k-i:-1:k+1-L-i);
    end
    
    hX = he_BSPAPSA(k,:)*Xap;
    
    h_bs=vec_div(he_BSPAPSA(k,:),div);
    h_n21=norm21(h_bs);%the norm 21s of the block || row vector
    lInf=max(abs(h_n21),delta);
    g=max(rou*lInf,abs(h_n21));
    gm=mean(g);
    g=g/gm;
            
    h_G=g(1)*ones(1,M);
    
    for j=2:div
        h_G=[h_G,g(j)*ones(1,M)];%build the block proportionate matrix
    end
    G=diag(h_G);
    
    d = h*Xap;
    evector = d-hX+vn(k:-1:k-L+1);
    e(k)=mean(evector);

    he_BSPAPSA(k+1,:)=he_BSPAPSA(k,:)+abs((abs(e(k))-abs(e(k-1)))/(abs(e(k))+abs(e(k-1))))*(mu*G*Xap*inv(Xap'*G*Xap+delta*I)*sign(evector'))';

end

for k=len/2:len-1
    for i=0:N-1
        Xap(i+1,:)=x(k-i:-1:k+1-L-i);
    end
    
    hX = he_BSPAPSA(k,:)*Xap;
    
    h_bs=vec_div(he_BSPAPSA(k,:),div);
    h_n21=norm21(h_bs);%the norm 21s of the block || row vector
    lInf=max(abs(h_n21),delta);
    g=max(rou*lInf,abs(h_n21));
    gm=mean(g);
    g=g/gm;
            
    h_G=g(1)*ones(1,M);
    
    for j=2:div
        h_G=[h_G,g(j)*ones(1,M)];%build the block proportionate matrix
    end
    G=diag(h_G);
    
    d = h1*Xap;
    evector = d-hX+vn(k:-1:k-L+1);
    e(k)=mean(evector);

    he_BSPAPSA(k+1,:)=he_BSPAPSA(k,:)+abs((abs(e(k))-abs(e(k-1)))/(abs(e(k))+abs(e(k-1))))*(mu*G*Xap*inv(Xap'*G*Xap+delta*I)*sign(evector'))';

end

he_BSPAPSA_final=he_BSPAPSA_final+he_BSPAPSA;
end
he_BSPAPSA_final=he_BSPAPSA_final/ROUND;

%%
% figure(1);
% plot(e(N:len).^2,'b');
% xlabel('迭代次数');
% ylabel('均方误差');
% title('APA')
% saveas(gcf,'C:\Work\毕设\codes\pic0409\APA_MSE.png')


figure(2)

Misalign_BSPAPA = 10*log10(sum((he_BSPAPA1_final-h).^2,2)/sum(h.^2));
Misalign_BSPAPA1 = 10*log10(sum((he_BSPAPA1_final-h1).^2,2)/sum(h1.^2));

Misalign_BSPAPSA = 10*log10(sum((he_BSPAPSA_final-h).^2,2)/sum(h.^2));
Misalign_BSPAPSA1 = 10*log10(sum((he_BSPAPSA_final-h1).^2,2)/sum(h1.^2));

plot([Misalign_BSPAPA(1:len/2-1);Misalign_BSPAPA1(len/2:len)],'k','linewidth',1.5)
hold on

plot([Misalign_BSPAPSA(1:len/2-1);Misalign_BSPAPSA1(len/2:len)],'r','linewidth',1.5)
hold on
legend('BS-PAPSA','EDTD-BS-PAPSA')
title('NMSE')
xlabel('迭代次数');
ylabel('归一化系数失调(dB)');