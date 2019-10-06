% clc;
% clear all;

clc;
clear;
tempp = 1;

tic

h_con = sqrt(3);
tf = 100; % 仿真时长
N = 50;  % 粒子个数

%系统矩阵设置
I = eye(3,3);                %生成单位阵
%I(3,3) = 0;

%噪声
q = [0 0 0]';                      %系统噪声均值
r = 1;                             %测量噪声均值
Q = [0.3,0,0; 0,0.3,0; 0,0,0.001]; %系统噪声方差矩阵
R = 0.3;                           %测量噪声方差矩阵
W_noise = Q*randn(3,N) + q;        %系统噪声
V_noise = R*randn(1,N) + r;        %测量噪声

x = zeros(3,tf); % 系统状态真实值 初始值0
P_pred_0 = [100,0,0; 0,100,0;0,0,100];%初始预测方差阵
x_0  = [0; 0; 0.2];                   %初始状态
x(:,1) = x_0;                    %真实状态初始值
y = zeros(1,tf); % 系统状态真实值 初始值0
y(1,1) = 0.1*x(1,1) + 0.3*x(2,1) +  sqrt(R) * randn + r;

P = cell(1,tf);        % 采样方差
P{1,1} = P_pred_0;     %初始估计方差阵
xhatPart = zeros(3,tf);%状态估计值

%记录运算时间以及误差的范数
FPF_SISO_ERROR_TIME = zeros(3,50);

xpart = zeros(3,N,tf);
for i = 1 : N
    xpart(:,i,1) = x(:,1) + sqrt(P{1,1}) * randn + q;%初始状态服从x=0均值，方差为sqrt(P)的高斯分布
end
% xArr = [x];
% yArr = [];
% xhatArr = [x];
% PArr = [P];
%xhatPartArr = [xhatPart]; %

%计算alpha阶次对应的GL定义系数 binomial coefficient
bino_fir = zeros(1,tf);       %微分阶次为0.7时GL定义下的系数
alpha = 1;
bino_fir(1,1) = 0.7;
for i = 2:1:tf
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);
end

bino_sec = zeros(1,tf);          %微分阶次为1.2时GL定义下的系数
alpha = 1.2;
bino_sec(1,1) = 1;
for i = 2:1:tf
    bino_sec(1,i) = (1-(alpha+1)/(i-1))*bino_sec(1,i-1);  
end

bino_thi = zeros(1,tf);          %微分阶次为1时GL定义下的系数
alpha = 1;
bino_thi(1,1) = 1;
for i = 2:1:tf
    bino_thi(1,i) = (1-(alpha+1)/(i-1))*bino_thi(1,i-1);  
end

%计算GL微分定义下系数矩阵
gamma = cell(1,tf);
temp_matrx = zeros(3,3);
for i = 1:1:tf 
    temp_matrx(1,1) = bino_fir(1,i);
    temp_matrx(2,2) = bino_sec(1,i);
    temp_matrx(3,3) = bino_thi(1,i);
    gamma{1,i} = temp_matrx;
end

%系统输入设置:前一半输入为1，后一半输入为
U_input = ones(1,tf);            
for i = N/2+1:1:tf
    U_input(:,i) = -2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diff_X_real    表示k时刻状态的微分
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff_X_real = [0;0;0];
xpartminus = zeros(3,N);

for k = 2 : tf
    
    diff_X_real = [x(2,k-1); ...
                  -0.1*x(1,k-1)-x(2,k-1)*x(3,k-1)+ ...
                  U_input(1,k-1); 0] + W_noise(:,k-1);
    %3*sin(2*x(1,k-1)) -x(1,k-1) + W_noise(1,k-1);
    rema = [0;0;0];
    for i = 2:1:k
        rema = rema + gamma{1,i}*x(1,k+1-i);
    end
    x(1,k) = diff_X_real - rema;
    %k时刻真实值
    y(1,k) = 0.1*x(1,k) + 0.3*x(2,k) + V_noise(1,k);  %k时刻观测值

 %% 采样N个粒子
 for i = 1 : N
     %采样获得N个粒子
     xpartminus(:,i) = [xpart(2,i,k-1); ...
                        -0.1*xpart(1,i,k-1)-xpart(2,i,k-1)*xpart(3,i,k-1)+ ...
                        U_input(1,k-1); 0] + sqrt(Q) * randn;
     %3*sin(2*xpart(i,k-1)) - xpart(i,k-1) + sqrt(Q) * randn;
     temp = 0;
         for j = 2 : 1 : k
            temp = temp + gamma{1,i}*xpart(:,i,k+1-j);
         end
     xpartminus(:,i) = xpartminus(:,i) - temp;
     ypart = xpartminus(:,i);      %每个粒子对应观测值
     vhat = y(1,k) - ypart;        %与真实观测之间的似然
     q(i) = (1 / sqrt(R) / sqrt(2*pi)) * exp(-vhat^2 / 2 / R);
     %每个粒子的似然即相似度
 end

 %%
%权值归一化
qsum = sum(q);
for i = 1 : N
    q(i) = q(i) / qsum; %归一化后的权值 q
end

%%
 %根据权值重新采样
  for i = 1 : N 
      u = rand;
      qtempsum = 0; 
      for j = 1 : N
          qtempsum = qtempsum + q(j); 
          if qtempsum >= u 
              xpart(i,k) = xpartminus(j);
              break;
          else
              xpart(i,k) = xpart(i,k-1);
          end
      end
  end
xhatPart(1,k) = mean(xpart(:,k));

%%
%最后的状态估计值即为N个粒子的平均值，这里经过重新采样后各个粒子的权值相同
% xArr = [xArr x];   
% yArr = [yArr y];  
% % xhatArr = [xhatArr xhat]; 
% PArr = [PArr P]; 
% xhatPartArr = [xhatPartArr xhatPart];

end

 FPF_SISO_TIME(1,tempp)    = toc
 FPF_ERROR_norm1(1,tempp)  = norm((x - xhatPart),1);
 FPF_ERROR_norm2(1,tempp)  = norm((x - xhatPart),2);
 tempp = tempp + 1;
 
%  FPF_SISO_ERROR_TIME(1,:)  = FPF_ERROR_norm1(1,:);
%  FPF_SISO_ERROR_TIME(2,:)  = FPF_ERROR_norm2(1,:);
%  FPF_SISO_ERROR_TIME(3,:)  = FPF_SISO_TIME(1,:);
%  FPF_ERROR_norm1_average  = sum(FPF_SISO_ERROR_TIME(1,:))/50
%  FPF_ERROR_norm2_average  = sum(FPF_SISO_ERROR_TIME(2,:))/50
%  FPF_SISO_TIME_average    = sum(FPF_SISO_ERROR_TIME(3,:))/50
%  save FPF_SISO_ERROE_TIME FPF_SISO_ERROR_TIME FPF_ERROR_norm1_average ...
%       FPF_ERROR_norm2_average FPF_SISO_TIME_average







