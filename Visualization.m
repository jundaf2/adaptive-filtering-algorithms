%parameters
N=3;L=4;len=500;
x=1*randn(1,len);
%%
d=[1 -1.5 0.7 0.1]; c=[1 0.5 0.2];  % 分子分母多项式系数
nd=length(d)-1 ;nc=length(c)-1;  %阶次
xik=zeros(nc,1);  %白噪声初值
ek=zeros(nd,1);
xi=randn(len,1);  %产生均值为0，方差为1的高斯白噪声序列

for k=1:len
   nColored(k)=-d(2:nd+1)*ek+c*[xi(k);xik];  %产生有色噪声
    %数据更新
    for i=nd:-1:2
       ek(i)=ek(i-1);
    end
   ek(1)=nColored(k);
    for i=nc:-1:2
       xik(i)=xik(i-1);
    end
   xik(1)=xi(k);
end
x=nColored;
vn = randn(1,len);
%%
varepsilon=0.01;
delta=0.01;
rou=0.01;
mu=0.01;
g(N)=0;
h=[0,2.5,0];

% X = zeros(1,N);
ROUND = 1;
I=eye(L);
Xap=zeros(N,L);

%% NLMS
he_NLMS = zeros(length(x),N);
% he_NLMS(:,1) = 2.5*ones(length(x),1);
% he_NLMS(:,3) = 2.5*ones(length(x),1);
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

%% APA
he_APA_final = zeros(length(x),N);
he_APA = zeros(length(x),N);
% he_APA(:,1) = 2.5*ones(length(x),1);
% he_APA(:,3) = 2.5*ones(length(x),1);
for round=1:ROUND
for k=N+L+1:len-1
    for i=0:N-1
        Xap(i+1,:)=x(k-i:-1:k+1-L-i);
    end
    hX = he_APA(k,:)*Xap;
    d = h*Xap;
    evector = d-hX+vn(k:-1:k-L+1);%randn(1,L);
    e(k)=mean(evector);

    he_APA(k+1,:)=he_APA(k,:)+(mu*Xap*inv(Xap'*Xap+delta*I)*evector')';
end
he_APA_final=he_APA_final+he_APA;
end
he_APA_final=he_APA_final/ROUND;
%% PAPA
F = @(x)min(2,400*x);

he_PAPA = zeros(length(x),N);
% he_PAPA(:,1) = 2.5*ones(length(x),1);
% he_PAPA(:,3) = 2.5*ones(length(x),1);
he_PAPA_final = zeros(length(x),N);
for round=1:ROUND
for k=N+L+1:len-1
    for i=0:N-1
        Xap(i+1,:)=x(k-i:-1:k+1-L-i);
    end
    hX = he_PAPA(k,:)*Xap;
    d = h*Xap;
    evector = d-hX+vn(k:-1:k-L+1);%randn(1,L);
    e(k)=mean(evector);
    lInf=max(F(abs(he_PAPA(k,:))),delta);
    g=max(rou*lInf,F(abs(he_PAPA(k,:))));
    gm=mean(g);
    g=g/gm;
    G=diag(g);
    he_PAPA(k+1,:)=he_PAPA(k,:)+(mu*G*Xap*inv(Xap'*G*Xap+delta*I)*evector')';
end
he_PAPA_final=he_PAPA_final+he_PAPA;
end
he_PAPA_final=he_PAPA_final/ROUND;
%%
% figure(1);
% plot(e(N:len).^2,'b');
% xlabel('迭代次数');
% ylabel('均方误差');
% title('APA')
% saveas(gcf,'C:\Work\毕设\codes\pic0409\APA_MSE.png')



d1=sum((he_NLMS_final-h).^2,2);
d2=sum((he_APA_final-h).^2,2);
d3=sum((he_PAPA_final-h).^2,2);
%%
figure(2)

hold on;

scatter3(h(1),h(2),h(3),10,'g','filled');
xlabel('w_0');
ylabel('w_1');
zlabel('w_2');
xlim([-0.5,3]);%对X轴设定显示范围 
ylim([-0.5,3]);%对Y轴设定显示范围 
zlim([-0.5,3]);%对Z轴设定显示范围 
i=1;
scatter3(he_NLMS_final(i,1),he_NLMS_final(i,2),he_NLMS_final(i,3),5,'k','filled')
hold on
scatter3(he_APA_final(i,1),he_APA_final(i,2),he_APA_final(i,3),5,'b','filled')
hold on
scatter3(he_PAPA_final(i,1),he_PAPA_final(i,2),he_PAPA_final(i,3),5,'r','filled')
hold on
h1=text(he_NLMS_final(i,1)+0.1,he_NLMS_final(i,2)+0.1,he_NLMS_final(i,3)-0.1,num2str(d1(i)),'color','k');
h2=text(he_APA_final(i,1)+0.1,he_APA_final(i,2)-0.1,he_APA_final(i,3)+0.1,num2str(d2(i)),'color','b');
h3=text(he_PAPA_final(i,1)-0.1,he_PAPA_final(i,2)+0.1,he_PAPA_final(i,3)+0.1,num2str(d3(i)),'color','r');
pause;
for i=2:len
    scatter3(he_NLMS_final(i,1),he_NLMS_final(i,2),he_NLMS_final(i,3),5,'k','filled')
    hold on
    scatter3(he_APA_final(i,1),he_APA_final(i,2),he_APA_final(i,3),5,'b','filled')
    hold on
    scatter3(he_PAPA_final(i,1),he_PAPA_final(i,2),he_PAPA_final(i,3),5,'r','filled')
    hold on
    grid on
    
    set(h1,'visible','off')
    set(h2,'visible','off')
    set(h3,'visible','off')
    
    h1=text(he_NLMS_final(i,1)+0.1,he_NLMS_final(i,2)+0.1,he_NLMS_final(i,3)-0.1,num2str(d1(i)),'color','k');
    h2=text(he_APA_final(i,1)+0.1,he_APA_final(i,2)-0.1,he_APA_final(i,3)+0.1,num2str(d2(i)),'color','b');
    h3=text(he_PAPA_final(i,1)-0.1,he_PAPA_final(i,2)+0.1,he_PAPA_final(i,3)+0.1,num2str(d3(i)),'color','r');
    
    
    
    legend('w_o','NLMS','APA','PAPA')
   drawnow; 
end
set(h1,'visible','on')
set(h2,'visible','on')
set(h3,'visible','on')

% saveas(gcf,'C:\Work\毕设\codes\pic0409\APA_Misalignment.png')





