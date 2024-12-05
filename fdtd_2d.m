function [u] = fdtd_2d(forward, srcx, srcy)
% GPR二维正演程序，采用CPML吸收边界的FDTD算法
% forward 正演所需的结构体
% src激励源位置（x，y）
% u接收信号
%%从结构体中获取正演参数
xdim = forward. xdim;
ydim = forward. ydim;
npml = forward. npml;
numit = forward. numit;
srcpulse = forward. srcpulse;
ds = forward.ds;
CA = forward.CA;
CB =forward.CB;
CQ= forward.CQ;
iEKx = forward. iEKx;
iHKx = forward. iHKx;
iEKy =forward. iEKy;
iHKy =forward. iHKy;
Bx_E =forward. Bx_E;
Ax_E =forward. Ax_E;
By_E =forward. By_E;
Ay_E = forward. Ay_E;
Bx_H =forward. Bx_H;
Ax_H = forward. Ax_H;
By_H =forward. By_H;
Ay_H= forward. Ay_H;
%%初始化%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u =zeros(xdim, ydim, numit);
PEzx = zeros(xdim, ydim);
PEzy = zeros(xdim, ydim);
PHx = zeros(xdim, ydim-1);
PHy = zeros(xdim -1, ydim);
Ez= zeros(xdim, ydim);
Hx = zeros(xdim, ydim-1);
Hy = zeros(xdim -1, ydim);
Ezdiffx = zeros( xdim -1, ydim);
Ezdiffy = zeros( xdim, ydim -1);
Hydiffx = zeros(xdim -1, ydim-1);% Hy在x方向的导数矩阵
Hxdiffy =zeros(xdim-1, ydim-1);% Hx在y方向的导数矩阵
%%时间循环%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1: numit

%%%%%%  HX  %%%%%%%%%
i=1: xdim;
j=1: ydim-1;

Ezdiffy =(Ez(i,j+1)- Ez(i,j))/ds;
Hx(i,j)= Hx(i, j)- CQ(i, j).* iHKy(i, j).* Ezdiffy(i, j);
% psi 上下
i=[1: npml, (xdim-npml +1):xdim];
j=1: ydim-1;
PHx(i,j) =By_H(i, j).* PHx(i, j)+ Ay_H(i, j).* Ezdiffy(i, j);
Hx(i,j)= Hx(i,j)- CQ(i,j).* PHx(i,j);

% 求psi（左右PML区域)
i=(npml +1): (xdim-npml); j=[1: npml, (ydim -npml +1): (ydim-1)];
PHx(i,j)= By_H(i, j).* PHx(i, j) + Ay_H(i,j).* Ezdiffy(i,j);
Hx(i,j)= Hx(i,j)- CQ(i,j).*PHx(i,j);
%%%%%%%%%%%%%%%%   HY %%%%%%%%%%%%%%%%5
i=1: (xdim-1); j=1: ydim;
Ezdiffx(i,j)=(Ez(i+1,j) - Ez(i, j))/ds;
Hy(i,j)= Hy(i,j)+CQ(i,j).* iHKx(i,j).* Ezdiffx(i, j);
%psi sx
i=[1: npml, (xdim -npml +1): (xdim-1)]; j=1: ydim;
PHy(i,j) = Bx_H(i, j).* PHy(i, j) + Ax_H(i, j).* Ezdiffx(i, j);
Hy(i,j) = Hy(i,j) + CQ(i, j).*PHy(i,j);
%%%%%%
i=(npml +1): (xdim-npml); j=[1: npml, (ydim -npml +1): ydim];
PHy(i,j)= Bx_H(i, j).* PHy(i, j) + Ax_H(i, j).* Ezdiffx(i, j);
Hy(i,j)= Hy(i,j)+ CQ(i,j).*PHy(i,j);
%%%%%%%%%
i=2: (xdim-1);j=2: (ydim-1);
Hydiffx(i,j)=(Hy(i, j) - Hy(i-1,j))/ds;
Hxdiffy(i,j)=(Hx(i,j)- Hx(i,j-1))/ds;
Ez(i,j)= CA(i,j).*Ez(i,j)+ CB(i,j).*( iEKx(i, j).* Hydiffx(i,j)-iEKy(i,j).* Hxdiffy(i,j));
%%%%%%%%%%%
i=[2: npml, (xdim -npml +1): (xdim-1)]; j=2: ydim-1;
PEzx(i, j)= Bx_E(i, j).* PEzx(i, j) + Ax_E(i, j).* Hydiffx(i, j);
PEzy(i, j)= By_E(i, j).* PEzy(i, j) + Ay_E(i, j).* Hxdiffy(i, j);
Ez(i,j)= Ez(i, j)+ CB(i, j).*(PEzx(i, j) - PEzy(i,j));
%%%%%%%%%
i=(npml +1): (xdim-npml); j=[2: npml,( ydim -npml +1): (ydim-1)];
PEzx(i,j)= Bx_E(i,j).* PEzx(i, j) + Ax_E(i, j).* Hydiffx(i, j);
PEzy(i,j)= By_E(i, j).* PEzy(i, j) + Ay_E(i, j).* Hxdiffy(i, j);
%把两个psi的值加入Ez中（左右PML区域)
Ez(i,j)= Ez(i,j)+ CB(i,j).*(PEzx(i, j) - PEzy(i, j));
Ez(srcx, srcy) = Ez(srcx, srcy) -CB(srcx, srcy) * srcpulse(n)/ds*2;  %加入激励源
u(:,:,n)=Ez;
end

