%
%   分数阶卡尔曼滤波器仿真复现
%   论文：     fractional order CDKF
%   目的：一阶分数阶中心差分卡尔曼滤波器算法测试
%         对系统噪声均值进行估计
%         函数实验:    D^{0.7} x_k = 3*sin(2*x_{k-1}) -x_{k-1} + w_k
%                              y_k = x_k + v_k
%   结果：
%
%   备注：分数阶中心差分卡尔曼滤波器算法
%         *********运算时间测试*********
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%  clc
%  clear
%  tempp = 1;

tic
 
%仿真步长
N = 100;
h_con = sqrt(2.0);

%状态测量初始化
X_real = zeros(1,N);         %真实状态
X_esti = zeros(1,N);         %状态最优估计值
P_esti = zeros(1,N);         %估计误差方差阵
Z_meas = zeros(1,N);         %实际观测值

%系统矩阵设置
I = eye(1,1);                %生成单位阵
 
%噪声
q = 0;                      %系统噪声均值
r = 0;                      %测量噪声均值
Q = 0.81;                   %系统噪声方差矩阵
R = 0.25;                   %测量噪声方差矩阵
W_noise = sqrt(Q)*randn(1,N) + q;  %系统噪声
V_noise = sqrt(R)*randn(1,N) + r;  %测量噪声

%初始值设置（初始矩阵不能为零）
P_pred_0 = 100;              %初始预测方差阵
P_esti(1,1) = P_pred_0;      %初始估计方差阵
x_0  = 0;                    %初始状态     
X_real(1,1) = x_0;           %真实状态初始值
X_esti(1,1) = 0;             %状态估计初始值
Z_meas(1,1) = V_noise(1,1);  %测量数据初始值

%FCDKF_SE = zeros(1,N);
FCDKF_SISO_ERROR_TIME = zeros(3,50);%一共做50次仿真运算

%计算alpha阶次对应的GL定义系数 binomial coefficient 
bino_fir = zeros(1,N);       %微分阶次为0.7时GL定义下的系数
alpha = 1;
bino_fir(1,1) = 0.7;
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
        diff_X_esti = 3*sin(2*X_esti(1,k-1)) - X_esti(1,k-1);
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
        X_pre = diff_X_esti - rema + q;         %一步状态预测
        %预测误差协方差矩阵:P_pred

            %对误差矩阵进行cholsky分解
            S_chol = chol(P_esti(1,k-1))';

            %计算余项
            rema_P = 0;
            if k>L+1
                for i = 2:1:L+2
                    rema_P = rema_P + bino_fir(1,i)*P_esti(1,k+1-i)*bino_fir(1,i)';
                end
            else
                for i = 2:1:k
                    rema_P = rema_P + bino_fir(1,i)*P_esti(1,k+1-i)*bino_fir(1,i)';
                end
            end

        %临时变量 temp_fun : 函数差值,函数为单变量
        temp_fun = 3*sin( 2*(X_esti(1,k-1)+h_con*S_chol) ) - (X_esti(1,k-1)+h_con*S_chol) - ...
                   3*sin( 2*(X_esti(1,k-1)-h_con*S_chol) ) + (X_esti(1,k-1)-h_con*S_chol);
        P_pred = 1/(4*h_con^2) * temp_fun^2 + rema_P + ...
                 1/h_con*0.5*temp_fun*S_chol'*(-bino_fir(1,2))' + ...
                 1/h_con*(-bino_fir(1,2))*S_chol*0.5*temp_fun' + Q;

        %测量值估计  Z_esti ---- Z_k|k-1
        Z_esti = X_pre + r;
        
        %测量预测误差协方差矩阵:P_zpred ---- P_z_k|k-1
        P_zpred = P_pred + R;

        %计算卡尔曼增益:Kk(2*1)
        Kk = P_pred/P_zpred;
        
        %状态更新
        X_esti(1,k) = X_pre + Kk*( Z_meas(1,k) - Z_esti );

        %估计方差矩阵更新
        P_esti(1,k) = P_pred - Kk*P_zpred*Kk';  
end

 FCDKF_SISO_TIME(1,tempp)    = toc
 FCDKF_ERROR_norm1(1,tempp)  = norm((X_real - X_esti),1)  
 FCDKF_ERROR_norm2(1,tempp)  = norm((X_real - X_esti),2)
 tempp = tempp + 1;

%  FCDKF_SISO_ERROR_TIME(1,:) = FCDKF_ERROR_norm1(1,:);
%  FCDKF_SISO_ERROR_TIME(2,:) = FCDKF_ERROR_norm2(1,:);
%  FCDKF_SISO_ERROR_TIME(3,:) = FCDKF_SISO_TIME(1,:);
%  FCDKF_ERROR_norm1_average  = sum(FCDKF_SISO_ERROR_TIME(1,:))/50
%  FCDKF_ERROR_norm2_average  = sum(FCDKF_SISO_ERROR_TIME(2,:))/50
%  FCDKF_SISO_TIME_average    = sum(FCDKF_SISO_ERROR_TIME(3,:))/50
%  save FCDKF_SISO_ERROE_TIME1 FCDKF_SISO_ERROR_TIME FCDKF_ERROR_norm1_average ...
%       FCDKF_ERROR_norm2_average FCDKF_SISO_TIME_average






