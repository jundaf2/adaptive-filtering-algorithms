%'PAPA','BS-PAPA' 之间的对比 。 所有程序 的 变量名称 不都完全和实际意义相同
N=512;L=4;len=200000; 
x=10*randn(1,len);

a1=-0.5;
for k=2:len
     x(k)=-a1*x(k-1)+x(k);
end 

varepsilon=0.01;
delta=0.01;
rou=0.01;
mu=0.01;
% h=load('randomGenerated_twoCluster_sparsEchopath.mat');
% h=h.w0;
% h=h';

h=zeros(1,N);
spst=12;
h(128:128-1+spst)=5*randn(1,spst);
h(256:256-1+spst)=5*randn(1,spst);

ROUND=1;
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
vn = randn(1,len);%2*iws(a,gama,ms,len);
figure (1);
plot(x);
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
    evector = d-hX+vn(k:-1:k-L+1);
    lInf=max(abs(he_PAPA(k,:)),delta);
    g=max(rou*lInf,abs(he_PAPA(k,:)));
    gm=mean(g);
    g=g/gm;
    G=diag(g);
    he_PAPA(k+1,:)=he_PAPA(k,:)+(mu*G*Xap*inv(Xap'*G*Xap+delta*I)*evector')';
end
he_PAPA_final=he_PAPA_final+he_PAPA;
end
he_PAPA_final=he_PAPA_final/ROUND;

%% BSPAPA 1
div=128; M=N/div;
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
    
    he_BSPAPA1(k+1,:)=he_BSPAPA1(k,:)+(mu*G*Xap*inv(Xap'*G*Xap+delta*I)*evector')';
end

he_BSPAPA1_final=he_BSPAPA1_final+he_BSPAPA1;
end
he_BSPAPA1_final=he_BSPAPA1_final/ROUND;

%%
figure(2)

Misalign_PAPA = 10*log10(sum((he_PAPA_final-h).^2,2)/sum(h.^2));

Misalign_BSPAPA1 = 10*log10(sum((he_BSPAPA1_final-h).^2,2)/sum(h.^2));

plot(Misalign_PAPA,'k','linewidth',2)
hold on
plot(Misalign_BSPAPA1,'b','linewidth',2)
hold on

legend('PAPA','BS-PAPA')
title('NMSE')
xlabel('迭代次数');
ylabel('归一化系数失调(dB)');