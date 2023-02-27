close all;
%% 参数
N = 100000;%输入信号长度
L = 128;%滤波器阶数
pi = 3.14;%圆周率
M = 8;%投影阶数
varepsilon = 0.01;%防止divide error
delta = 0.01;%防止所有系数都为0时算法冻结（迭代增量为0）
% rou = 5/L;%防止某个系数过小时该系数不更新（成比例分配的比例为0）
theta = 0;%系统稀疏度
alpha = 1-M/(4*L);%平滑因子
Nw = 16;%用于噪声估计的估计窗
C = 1.483*(1+5/(Nw-1));%有限样值矫正因子

%% 仿射投影矩阵
Xap = zeros(M,L);

%% 定义向量
e = zeros(N,M);
Ae = zeros(L,Nw);
Be = zeros(L,Nw);
beta = zeros(1,N);
sigma_x = zeros(1,N);
sigma_e = zeros(1,N);
sigma_v = zeros(1,N);
r = zeros(1,Nw);
%% 信道模型
k = 1:L;
h = 1./k.*exp(-(k-L/2).^2/8000).*sin(pi*(k-1)/6);
figure(1)
plot(h);title('非稀疏信道');hold on;

he = zeros(N,L);% 信道估计矩阵，记录迭代过程中的变化

%% 输入信号为色噪声
num = [1];
den = [1,-0.9];
x = filter(num,den,randn(1,N));
figure(2)
plot(x);title('有色信号');
x=10*randn(1,N);

%% alpha-stable clutter 噪声时间序列
a = 1.5;gama = 1.5;ms = 1;
vn = iws(a,gama,ms,N);
figure(3)
plot(vn);
title('噪声v');
%% 变步长
mu = ones(1,N);

%% 比例系数矩阵
g = ones(1,L);
t=1;

%% 迭代
mu0 = norm(h,inf)/20;

sigma_x(L+M-1) = var(x(1:L+M));
for k = L+M:N
    for i = 1:L
        Xap(:,i) = x(k+1-i:-1:k+2-M-i)';
    end
    g = (1-theta)/(2*L)+(1+theta)*abs(he(k,:))/(2*norm(he(k,:),1)+varepsilon);
    G = diag(g);
    e(k,:) = (Xap*h'-Xap*he(k,:)')'+vn(k:-1:k-M+1);
    if mod(k-L-M,Nw)==0
        t=t+1;
        Ae = e(k:-1:k-Nw+1,:)';
        Be = e(k:-1:k-Nw+1,:).^2';
        r = alpha*r+C*(1-alpha)*x(k)*median(Ae);
        sigma_x(t) = alpha*sigma_x(t-1)+(1-alpha)*mu(t-1)^2*0.5;
        sigma_e(t) = alpha*sigma_e(t-1)+C*(1-alpha)*median(median(Be))*0.5;
        sigma_v(t) = sigma_e(t)+r*r'/sigma_x(t)*0.5;
        beta(t) = (norm(e(t,:),1)-M*sqrt(2/pi*sigma_v(t)))/sqrt(sign(e(t,:))*(Xap*Xap')*sign(e(t,:)'));
        SignOfBeta = [mu(t-1),0,beta(t)];%避免if else
        beta(t) = SignOfBeta(sign(beta(t))+2);%把负数映射为下标1，把正数映射为下标3
        mu(t) = alpha*mu(t-1)+(1-alpha)*min(beta(t),mu(t-1));
    end
    he(k+1,:) = he(k,:)+mu0*mu(t)*sign(e(k,:))*Xap*G*inv(sqrt(sign(e(k,:))*Xap*G*Xap'*sign(e(k,:)')));
%      figure(1);
%      fig=plot(he(k,:),'r');title('信道估计');drawnow;
%      set(fig,'visible','off');
%      figure(3)
%      plot([vn(1:k),zeros(1,N-k)]);drawnow;
end

%%
figure(4);
plot(mean(e(1:N,:),2).^2,'b');
xlabel('迭代次数');
ylabel('均方误差');
title('VSS_IPAPSA')
% saveas(gcf,'C:\Work\毕设\codes\pic0409\VSS_IPAPSA_MSE.png')

figure(5)
Misalign_VSS_IPAPSA = 10*log10(sum((he-h).^2,2)/sum(h.^2));
plot(Misalign_VSS_IPAPSA)
xlabel('迭代次数');
ylabel('系数失调(dB)');
title('VSS_IPAPSA');
% saveas(gcf,'C:\Work\毕设\codes\pic0409\VSS_IPAPSA_Misalignment.png')
% 
% figure(6)
% plot(he(N,:))
% 
