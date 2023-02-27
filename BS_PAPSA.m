%parameters
N=512;L=8;len=20000; div=32; M=N/div;
x=10*randn(1,len);
varepsilon=0.01;
delta=0.01;
rou=0.01;
mu=1;
h=load('randomGenerated_twoCluster_sparsEchopath.mat');
h=h.w0;
h=h';
he = zeros(length(x),N);
% X = zeros(1,N);
v = randn(1,length(x));
var_x = zeros(1,length(x));
e = zeros(1,length(x));
I=eye(L);
Xap=zeros(N,L);
h_bs=zeros(div,N/div);
IM=eye(M);
% F is a real-valued function to map the current coe?cient estimate into 
% a certain value of the proportionate step-size parameter
F = @(x)min(2,400*x);


for k=N+L+1:len
    for i=0:N-1
        Xap(i+1,:)=x(k-i:-1:k+1-L-i);
    end
    
    hX = he(k,:)*Xap;
    
    h_bs=vec_div(he(k,:),div);
    h_n21=norm21(h_bs);%the norm 21s of the block || row vector
    lInf=max(F(abs(h_n21)),delta);
    g=max(rou*lInf,F(abs(h_n21)));
    gm=mean(g);
    g=g/gm;
            
    h_G=g(1)*ones(1,M);
    
    for j=2:div
        h_G=[h_G,g(j)*ones(1,M)];%build the block proportionate matrix
    end
    G=diag(h_G);
    
    d = h*Xap;
    evector = d-hX+randn(1,L);
    e(k)=mean(evector);
    
    he(k+1,:)=he(k,:)+(mu*G*Xap*inv(Xap'*G*Xap+delta*I)*sign(evector)')';
end
%%
figure(1);
plot(e(N:len).^2,'b');
xlabel('迭代次数');
ylabel('均方误差');
title('PAPA')
% saveas(gcf,'C:\Work\毕设\codes\pic0409\PAPA_MSE.png')


figure(2)
Misalign_PAPA = 10*log10(sum((he-h).^2,2)/sum(h.^2));
plot(Misalign_PAPA)
xlabel('迭代次数');
ylabel('系数失调(dB)');
title('PAPA')
% saveas(gcf,'C:\Work\毕设\codes\pic0409\PAPA_Misalignment.png')





