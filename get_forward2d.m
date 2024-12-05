function [forward] = get_forward2d(ep, sig, freq, ds,t, npml)
%光速c0
%真空中的介电常数eps0
%真空中的磁导率mu0
% PML层中的指数系数4
%生成GPR二维正演程序CPML吸收边界的FDTD算法所需系数的结构体
%输入：
% ep介电常数矩阵
% sig电导率矩阵 
% freq 主频
% ds空间间隔
% t时间向量
% npml pml 层数
%输出
% forward 正演所需系数结构体%%%%%%%%%%%%%%%%
c0= 2.99792458e8;
eps0=8.854187818e-12;
mu0=4*pi*1e-7;
m=4;
[xdim, ydim] =size(ep);
numit = length(t);
%模拟区域物性参数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps1 =ep;
sigma =sig;
mul = ones(xdim, ydim);
v=c0./sqrt(eps1(:,:).*mul(:,:));
c_max =max(max(v));
dt =t(2)-t(1);
%%%%规范时间步长与网格大小

dt_limit=ds/(sqrt(2)*c_max);
if dt > dt_limit
disp('dt过大')
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%布莱克曼一哈里斯脉冲
srcpulse = blackharrispulse(freq,t);
%初始化
EKx = ones(1, xdim);
HKx = ones(1,xdim);
iEKx = ones(1, xdim);
iHKx = ones(1, xdim);
Esigx = zeros(1, xdim);
Hsigx = zeros(1, xdim);
Ealphax = zeros(1, xdim);
Halphax = zeros(1, xdim);
EKy = ones(1, ydim);
HKy = ones(1, ydim);
iEKy = ones(1, ydim);
iHKy = ones(1, ydim);
Esigy = zeros(1, ydim);
Hsigy = zeros(1, ydim);
Ealphay = zeros(1, ydim);
Halphay = zeros(1, ydim);
%%%%%%%%%%%%-------------------------------------------%%%%%%%%%%%%%%%%%%%%%%%
%CPML

sigmax =(m +1)/(150* sqrt(eps1(1,1))* pi* ds);% PML 层的最大电导率
Kmax = 5;% kappa 最大值，kappa值的引入是为了改善PML对表面波的吸收特性
alphamax =0.008;% 0-0.05，alpha 值是为了改善PML对低频分量的
%左侧PML
i=1:npml;
Esigx(1, i) = sigmax *(((npml +1)-i)/npml).^m;
EKx(1,i)=1+(Kmax-1)*(((npml+1)-i)/npml).^m;
Ealphax(1, i)=alphamax*((i)/(npml));

Hsigx(1, i) =sigmax *((i-(npml+0.5))/npml).^m; %(10.5-i)*ds
HKx(1,i)=1+(Kmax-1)*((i-(npml+0.5))/npml).^m;
Halphax(1, i) =alphamax *((i +0.5)/(npml));
%右侧 PML
i=xdim -npml +1: xdim;
Esigx(1,i)= sigmax *((i -(xdim -npml))/npml).^m;
EKx(1,i)=1 +(Kmax -1)*((i-(xdim -npml))/npml).^m;
Ealphax(1, i) =alphamax *((xdim +1 -i)/(npml));

Hsigx(1, i-1) = sigmax * ((i-1 -(xdim -npml -0.5))/npml).^m; %(49.5-i)* ds
HKx(1, i-1) =1+(Kmax - 1)*((i-1-(xdim-npml-0.5))/npml).^m;
Halphax(1, i-1) =alphamax *((xdim +0.5-i+1)/(npml));
%上方 PML
j=1: npml;
Esigy(1, j) =sigmax *(((npml +1)-j)/npml).^m;
EKy(1,j)=1+(Kmax-1)*(((npml+1)-j)/npml).^m;
Ealphay(1, j) =alphamax*((j)/(npml));

Hsigy(1, j) = sigmax *(((npml +0.5)-j)/npml).^m;
HKy(1,j)=1+(Kmax-1)*(((npml +0.5)-j)/npml).^m;
Halphay(1, j)=alphamax *((j+0.5)/(npml));
%下方 PML
j =ydim -npml +1: ydim;
Esigy(1, j)=sigmax*((j-(ydim-npml))/npml).^m;
EKy(1,j)=1+(Kmax-1)*((j-(ydim-npml))/npml).^m;
Ealphay(1,j) =alphamax *((ydim +1-j)/(npml));

Hsigy(1,j-1) =sigmax *((j-1-(ydim-npml-0.5))/npml).^m;
HKy(1,j-1)=1 +(Kmax-1)*((j-1-(ydim-npml-0.5))/npml).^m;
Halphay(1, j-1) =alphamax*((ydim +0.5-j+1)/(npml));
%倒数
iEKx =1./EKx;
iHKx =1./HKx;
iEKy =1./EKy;
iHKy =1./HKy;
%分配Psi系数内存
Bx_E= zeros(1, xdim);
Ax_E= zeros(1, xdim);
By_E= zeros(1, ydim);
Ay_E= zeros(1, ydim);
By_H= zeros(1, ydim);
Ay_H= zeros(1, ydim);
Bx_H= zeros(1, xdim);
Ax_H= zeros(1,xdim);
for i=1: xdim
if (i <= npml)||(i >= xdim -npml+1)
%求Psi的第一个和第二个系数
Bx_E(1,i) = exp( -(dt/eps0)*((Esigx(1,i)./EKx(1,i))+Ealphax(1,i)));
Ax_E(1,i)=(Esigx(1,i)*( Bx_E(1,i) -1))/...
    (EKx(1,i) *( Esigx(1, i) + Ealphax(1,i) * EKx(1,i)));
end
if (i <= npml)||(i>=xdim -npml)
Bx_H(1,i) = exp(-(dt/eps0)*(Hsigx(1,i)./HKx(1,i)+Halphax(1,i)));
Ax_H(1,i)=(Hsigx(1,i)*( Bx_H(1,i) - 1))/...
    (HKx(1,i)*( Hsigx(1, i) + Halphax(1, i) * HKx(1,i)));
end
end
for j=1: ydim
if (j <= npml)||(j>= ydim-npml +1)
By_E(1,j) = exp( -(dt/eps0) *(Esigy(1, j)./EKy(1, j)+ Ealphay(1,j)));
Ay_E(1,j)=(Esigy(1,j)*( By_E(1,j) - 1))/...
    (EKy(1,j)*( Esigy(1,j) + Ealphay(1,j) * EKy(1,j)));
end
if (j <= npml)||(j>= ydim - npml)
By_H(1,j)= exp( -(dt/eps0)*(Hsigy(1,j)./HKy(1,j)+Halphay(1,j)));
Ay_H(1,j)=(Hsigy(1,j)*( By_H(1,j)- 1))/...
(HKy(1,j)*(Hsigy(1, j) + Halphay(1, j)* HKy(1,j)));
end
end
%%--------------------------------------------------------------------------------%%
CA =zeros(xdim, ydim); CB =zeros( xdim, ydim); CQ=zeros( xdim, ydim);
i=1: xdim;j=1: ydim;
CA(i, j)=(2*eps1(i, j) * eps0 -sigma(i, j)* dt)./(2*eps1(i, j)*eps0+ sigma(i, j)* dt); %CA
CB(i, j) =2* dt./(2 *eps1(i, j) * eps0 + sigma(i, j)* dt); %CB
CQ(i, j) = dt./(mu0 * mul(i, j));%CQ
%%--------------------------------------------------------------------------------%%
% repmat(将一维矩阵变成二维矩阵)
iEKx =repmat(iEKx', 1, ydim);
iHKx =repmat(iHKx', 1, ydim);
iEKy =repmat(iEKy, xdim, 1);
iHKy = repmat(iHKy,xdim,1);
% Psi系数矩阵（将一维矩阵变为二维矩阵）
Bx_E = repmat( Bx_E', 1, ydim);
Ax_E = repmat( Ax_E', 1, ydim);
By_E = repmat( By_E, xdim, 1);
Ay_E = repmat(Ay_E, xdim, 1);
Bx_H = repmat( Bx_H', 1, ydim);
Ax_H =repmat( Ax_H', 1, ydim);
By_H =repmat( By_H, xdim, 1);
Ay_H = repmat( Ay_H, xdim, 1);
%%生成正演所需结构体
forward. xdim = xdim;
forward. ydim = ydim;
forward.npml = npml;
forward. numit= numit ;
forward. srcpulse = srcpulse;
forward. ds= ds;
forward. CA = CA;
forward.CB = CB ;
forward. CQ = CQ;
forward. iEKx = iEKx;
forward. iHKx = iHKx;
forward. iEKy = iEKy;
forward. iHKy = iHKy;
forward. Bx_E = Bx_E;
forward. Ax_E = Ax_E;
forward. By_E = By_E;
forward. Ay_E= Ay_E;
forward. Bx_H = Bx_H;
forward. Ax_H = Ax_H;
forward. By_H = By_H;
forward. Ay_H = Ay_H;
