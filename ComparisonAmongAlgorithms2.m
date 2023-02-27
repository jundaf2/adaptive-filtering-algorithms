% 'PAPA','BS-PAPA, N=L=512','BS-PAPA, N=L/128=4' 不同块大小的BS-pAPA与 PAPA
N=512;L=2;len=5000; 
x=10*randn(1,len);
varepsilon=0.01;
delta=0.01;
rou=0.01;
mu=1;
h=load('randomGenerated_twoCluster_sparsEchopath.mat');
h=h.w0;
h=h';

h=zeros(1,length(h));
spst=100;
h(128:128-1+spst)=10;
h(290:290-1+spst)=10;

ROUND=1;
% X = zeros(1,N);
v = randn(1,length(x));
var_x = zeros(1,length(x));
I=eye(L);
Xap=zeros(N,L);

% F is a real-valued function to map the current coe?cient estimate into 
% a certain value of the proportionate step-size parameter
F = @(x)min(2,400*x);
%% PAPA
he_PAPA = zeros(length(x),N);
he_PAPA_final = zeros(length(x),N);
for round=1:ROUND
for k=N+L+1:len-1
    for i=0:N-1
        Xap(i+1,:)=x(k-i:-1:k+1-L-i);
    end
    hX = he_PAPA(k,:)*Xap;
    d = h*Xap;
    evector = d-hX+randn(1,L);
    lInf=max(F(abs(he_PAPA(k,:))),delta);
    g=max(rou*lInf,F(abs(he_PAPA(k,:))));
    gm=sum(g);
    g=g/gm;
    G=diag(g);
    he_PAPA(k+1,:)=he_PAPA(k,:)+(mu*G*Xap*inv(Xap'*G*Xap+delta*I)*evector')';
end
he_PAPA_final=he_PAPA_final+he_PAPA;
end
he_PAPA_final=he_PAPA_final/ROUND;

%% BSPAPA 1
div=256; M=N/div;
he_BSPAPA1 = zeros(length(x),N);
h_bs=zeros(div,M);
IM=eye(M);


he_BSPAPA1_final = zeros(length(x),N);
for round=1:ROUND
for k=N+L+1:len-1
    for i=0:N-1
        Xap(i+1,:)=x(k-i:-1:k+1-L-i);
    end
    
    hX = he_BSPAPA1(k,:)*Xap;
    
    h_bs=vec_div(he_BSPAPA1(k,:),div);
    h_n21=norm21(h_bs);%the norm 21s of the block || row vector
    lInf=max(F(abs(he_BSPAPA1(k,:))),delta);
    g=max(rou*lInf,F(abs(he_BSPAPA1(k,:))));
    gm=sum(g);
    g=g/gm;
            
    h_G=g(1)*ones(1,M);
    
    for j=2:div
        h_G=[h_G,g(j)*ones(1,M)];%build the block proportionate matrix
    end
    G=diag(h_G);
    
    d = h*Xap;
    evector = d-hX+randn(1,L);
    
    he_BSPAPA1(k+1,:)=he_BSPAPA1(k,:)+(mu*G*Xap*inv(Xap'*G*Xap+delta*I)*evector')';
end

he_BSPAPA1_final=he_BSPAPA1_final+he_BSPAPA1;
end
he_BSPAPA1_final=he_BSPAPA1_final/ROUND;


%% BSPAPA 2
div=4; M=N/div;
he_BSPAPA2 = zeros(length(x),N);
h_bs=zeros(div,M);
IM=eye(M);


he_BSPAPA2_final = zeros(length(x),N);
for round=1:ROUND
for k=N+L+1:len-1
    for i=0:N-1
        Xap(i+1,:)=x(k-i:-1:k+1-L-i);
    end
    
    hX = he_BSPAPA2(k,:)*Xap;
    
    h_bs=vec_div(he_BSPAPA2(k,:),div);
    h_n21=norm21(h_bs);%the norm 21s of the block || row vector
    lInf=max(F(abs(he_BSPAPA2(k,:))),delta);
    g=max(rou*lInf,F(abs(he_BSPAPA2(k,:))));
    gm=sum(g);
    g=g/gm;
            
    h_G=g(1)*ones(1,M);
    
    for j=2:div
        h_G=[h_G,g(j)*ones(1,M)];%build the block proportionate matrix
    end
    G=diag(h_G);
    
    d = h*Xap;
    evector = d-hX+randn(1,L);
    
    he_BSPAPA2(k+1,:)=he_BSPAPA2(k,:)+(mu*G*Xap*inv(Xap'*G*Xap+delta*I)*evector')';
end

he_BSPAPA2_final=he_BSPAPA2_final+he_BSPAPA2;
end
he_BSPAPA2_final=he_BSPAPA2_final/ROUND;

%%
% figure(1);
% plot(e(N:len).^2,'b');
% xlabel('迭代次数');
% ylabel('均方误差');
% title('APA')
% saveas(gcf,'C:\Work\毕设\codes\pic0409\APA_MSE.png')


figure(2)

Misalign_PAPA = 10*log10(sum((he_PAPA_final-h).^2,2)/sum(h.^2));

Misalign_BSPAPA1 = 10*log10(sum((he_BSPAPA1_final-h).^2,2)/sum(h.^2));

Misalign_BSPAPA2 = 10*log10(sum((he_BSPAPA2_final-h).^2,2)/sum(h.^2));
plot(Misalign_PAPA,'k','linewidth',1.5)
hold on
plot(Misalign_BSPAPA1,'b','linewidth',1.5)
hold on
plot(Misalign_BSPAPA2,'r','linewidth',1.5)
hold on
legend('PAPA','BS-PAPA, N=L=512','BS-PAPA, N=L/128=4')
title('NMSE')
xlabel('迭代次数');
ylabel('归一化系数失调(dB)');