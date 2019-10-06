%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   分数阶卡尔曼滤波器仿真复现
%   论文：     fractional order CDKF
%   目的：一阶分数阶中心差分卡尔曼滤波器算法测试
%         对系统噪声均值进行估计
%         函数实验:    D^{0.7} x_k = 3*sin(2*x_{k-1}) -x_{k-1} + w_k
%                              y_k = x_k + v_k
%   结果：较好的对状态进行估计，常值系统噪声均值收敛
%
%   备注：分数阶扩展卡尔曼滤波器的算法测试
%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  clc
%  clear
%  tempp = 1;

 tic

%仿真步长
N = 50;

%状态测量初始化
X_real = zeros(1,N);         %真实状态
X_esti = zeros(1,N);         %状态最优估计值
P_esti = zeros(1,N);         %估计误差方差阵
Z_meas = zeros(1,N);         %实际观测值

%系统矩阵设置
% A = [0,1; -0.1,-0.2];      %系统矩阵
% B = [0; 1];                %
% C = [0.1,0.3];             %
I = eye(1,1);                %生成单位阵
%I(3,3) = 0;

%噪声
q = 0;                     %系统噪声均值
r = 0;                     %测量噪声均值
Q = 0.81;                    %系统噪声方差矩阵
R = 0.25;                    %测量噪声方差矩阵
W_noise = sqrt(Q)*randn(1,N) + q;  %系统噪声
V_noise = sqrt(R)*randn(1,N) + r;  %测量噪声

%初始值设置
P_pred_0 = 100;              %初始预测方差阵
P_esti(1,1) = P_pred_0;      %初始估计方差阵
x_0  = 0;                    %初始状态
X_real(1,1) = x_0;           %真实状态初始值
X_esti(1,1) = 0;             %状态估计初始值
Z_meas(1,1) = V_noise(1,1);  %测量数据初始值

FEKF_SISO_ERROR_TIME = zeros(3,50);%一共做50次仿真运算

%计算alpha阶次对应的GL定义系数 binomial coefficient 
bino_fir = zeros(1,N);       %微分阶次为0.7时GL定义下的系数
alpha = 0.7;
bino_fir(1,1) = 1;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);  
end


%GL定义下短记忆原理的长度
L = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diff_X_real    表示k时刻状态的微分
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff_X_real = 0;


for k=2:1:N
    %计算实际状态
    diff_X_real = 3*sin(2*X_real(1,k-1)) -X_real(1,k-1) + W_noise(1,k-1);
    rema = 0;
    for i = 2:1:k
        rema = rema + bino_fir(1,i)*X_real(1,k+1-i);
    end
    X_real(1,k) = diff_X_real - rema;
    %实际观测值
    Z_meas(1,k) = X_real(1,k) + V_noise(1,k);
  %卡尔曼滤波
        %状态预测:X_pre
        diff_X_esti = 3*sin(2*X_esti(1,k-1)) - X_esti(1,k-1) + q;
            %计算余项
            rema = 0;
            if k>L
                for i = 2:1:L+1
                   rema = rema + bino_fir(1,i)*X_esti(1,k+1-i);
                end
            else
                for i = 2:1:k
                    rema = rema + bino_fir(1,i)*X_esti(1,k+1-i);
                end
            end
        X_pre = diff_X_esti - rema;     %一步状态预测
        %预测误差协方差矩阵:P_pred
            %计算余项
            rema_P = 0;
            if k>L+1
                for i = 3:1:L+2
                    rema_P = rema_P + bino_fir(1,i)*P_esti(1,k+1-i)*bino_fir(1,i)';
                end
            else
                for i = 3:1:k
                    rema_P = rema_P + bino_fir(1,i)*P_esti(1,k+1-i)*bino_fir(1,i)';
                end
            end
            F = 6*cos(2*X_esti(1,k-1)) - 1;
        P_pred = (F-bino_fir(1,2))*P_esti(1,k-1)*(F-bino_fir(1,2))'+ Q + rema_P;
        
        %计算卡尔曼增益:Kk(2*1)
        H = 1;
        Kk = P_pred*H'/(H*P_pred*H' + R);
        
        %状态更新
        X_esti(1,k) = X_pre + Kk*( Z_meas(1,k) - X_pre -r);
        
        %估计方差矩阵更新
        P_esti(1,k) = (I-Kk*H)*P_pred;

end

 FEKF_SISO_TIME(1,tempp)    = toc
 FEKF_ERROR_norm1(1,tempp)  = norm((X_real - X_esti),1)
 FEKF_ERROR_norm2(1,tempp)  = norm((X_real - X_esti),2)
 tempp = tempp + 1;
 
%  FEKF_SISO_ERROR_TIME(1,:) = FEKF_ERROR_norm1(1,:);
%  FEKF_SISO_ERROR_TIME(2,:) = FEKF_ERROR_norm2(1,:);
%  FEKF_SISO_ERROR_TIME(3,:) = FEKF_SISO_TIME(1,:);
%  FEKF_ERROR_norm1_average  = sum(FEKF_SISO_ERROR_TIME(1,:))/50
%  FEKF_ERROR_norm2_average  = sum(FEKF_SISO_ERROR_TIME(2,:))/50
%  FEKF_SISO_TIME_average    = sum(FEKF_SISO_ERROR_TIME(3,:))/50
%  save FEKF_SISO_ERROE_TIME1 FEKF_SISO_ERROR_TIME FEKF_ERROR_norm1_average ...
%       FEKF_ERROR_norm2_average FEKF_SISO_TIME_average










