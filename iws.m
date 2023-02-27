%% Generation of complex isotropic SaS random variables
% 【1】Simulation of Dependent Samples of Symmetric Alpha-Stable Clutter
% 【2】On the Chambers-Mallows-Stuck method for simulating skewed stable random variables
% 【3】The robust covariation-based MUSIC (ROC-MUSIC) algorithm for bearing estimation in impulsive noise environments
%% 此函数参考：文献【3】
function x=iws(alpha,gama,ms,ns)
% 0<alpha<=2:复SaS特征指数
% sigmaG:复Gaussian方差
%% 产生文献【1】式(4)：非相关正Alpha分布序列
p1=(2-alpha)/alpha;
rand('state',sum(100*clock));
vv=pi*(rand	(ms,ns));% 产生(0,pi)的均匀分布随机变量
rand('state',7*sum(100*clock)+3);
ww=-log(1-rand(ms,ns));% 产生均值为1的指数分布随机变量
yAw=sin(alpha*vv/2)./(sin(vv)).^(2/alpha).*(sin((1-alpha/2)*vv)./ww).^p1;% 产生nonnegative innovations η(k)
%% 产生文献【1】式(3)：复Gaussian分布序列
sigmaG=4;
sigmaG_dB=10*log10(sigmaG);
yGw=wgn(ms,ns,sigmaG_dB,'real');
xasigmac=std(yGw);
xamuc=mean(yGw);
yG=2*(yGw-xamuc)/xasigmac;
%% 产生非相关复SαS随机序列
x=gama^(1/alpha)*sqrt(yAw).*yG;
%% 估计非相关复SaS序列参数
sigmaA=sqrt(sigmaG)/2;
gama=sigmaA^alpha;



