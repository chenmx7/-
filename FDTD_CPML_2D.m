%FDTD(有限差分) 2D CPML
close all;
clear;
clc;
%%------设置模拟区域大小---------------------------%%
npml=10;  %PML边界层数
xmodel =200; %网格大小
ymodel =200;
xdim = xmodel +2 * npml;%x方向（纵向)总长
ydim = ymodel +2 * npml;%y方向（横向）总长
%% ——-设置天线参数以及空间时间步长-----------------%%
freq = 9e8;% 主频
ds =0.005; %网格步长
dt =1e-11; %时间步长
n_timestep= 1000; %1000个单位时
x=(1:1:xmodel) * ds;
y=(1:1:ymodel) * ds;
t=(0: 1: n_timestep-1) *dt; % 时间向量
%%----------信号源位置----------------------------%%
xsite = npml+1;
ysite = npml+1:5:ymodel+npml;
%% — ——--------模拟区域物性参数一—------------------%%

eps1 =6* ones( xdim, ydim);%相对介电常数
mul = ones(xdim, ydim); %相对磁导率

sigma =0.008 * ones(xdim,ydim); % 电导率

%%-------------------PML边界检验图--------------------%%
figure(1)
[forward0] = get_forward2d(eps1, sigma, freq, ds, t, npml); %获取正演结构体
[totalrecord] =fdtd_2d(forward0, xdim/2, ydim/2);
pic_num = 1;
for i=1:n_timestep
    imagesc(totalrecord(:,:,i));
    axis square;
    title(['CPML吸收边界条件','t = ' num2str(dt*i*1e9),'ns']);
    colorbar
    drawnow;
    F=getframe(gcf);
    I=frame2im(F);
    [I,map]=rgb2ind(I,256);
    if pic_num == 1
        imwrite(I,map,'CPML.test.gif','gif', 'Loopcount',inf,'DelayTime',0.01);
    else
        imwrite(I,map,'CPML.test.gif','gif','WriteMode','append','DelayTime',0.01);
    end
    pic_num = pic_num + 1;
end 

%%%%--------------------空洞模型-------------------------%%%
eps1(1:50,:)=16;
for i=1:xdim
    for j=1:ydim
        if((i-110)^2+(j-110)^2<=1000)
            eps1(i,j)=1;
            sigma(i,j)=0;
        end
    end
end
%%---------------------相对介电常数模型图-------------------%%

figure(2);
imagesc(x, y, eps1 (npml + 1: npml + xmodel, npml + 1: npml + ymodel));
ylabel('深度(m)');
axis image;
xlabel('水平位置(m)');
colorbar;
clim([1, 12])
title('相对介电常数模型图');

%%--------------自激自收(一点1000个)------------------%%
record=ones(n_timestep,length(ysite));
[forward] = get_forward2d(eps1, sigma, freq, ds, t, npml); %获取正演结构体
%%%-----------------道循环------------------%%
for i=1:length(ysite)

    [u] =fdtd_2d(forward, xsite, ysite);

record(:,i)=shiftdim(u(xsite,ysite(i),:));
disp(num2str(i));
end

toc
%%-------------------各点数据图---------------------%%
figure(3);
imagesc(y(1:5: end), t*1e9, record, [ -40, 40]);
xlabel('水平位置(m)');
ylabel('时间（ns）');
colorbar;
colormap jet
title('自激自收正演结果图');