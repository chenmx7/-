function srcpulse = blackharrispulse(freq, t)
%通过给定的时间向量和中心频率 freq计算 Blackman-Harris 窗函数的导
%输出： srcpulse=布菜克曼一哈里斯脉冲
%输入：freq=中心频率（Hz)
%t =时间向量（s）
%布莱克曼-哈里斯窗函数的系数
a= [0.35322222 -0.488 0.145 -0.010222222];
%改变T=1.14/freq使得频谱最大值接近中心频率
T= 1.14/freq;
f= zeros( size(t));
for n=0:3
f= f + a(n+1)*cos(2*n*pi*t./T);
end
f(t>=T) = 0;
%计算窗函数的一阶导数，并进行归一化
srcpulse= [f(2: end),0] -f(1: end);
srcpulse = srcpulse./max(abs( srcpulse));
