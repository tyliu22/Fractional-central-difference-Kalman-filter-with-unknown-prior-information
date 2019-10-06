%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   自适应分数阶中心差分卡尔曼滤波器估计性能分析
%   论文：     fractional order CDKF
%   目的：每次生成100次Q和r的估计数据，并保存至 estiqr.mat
%         函数实验:    D^{0.7} x_k = 3*sin(2*x_{k-1}) -x_{k-1} + w_k
%                              y_k = x_k + v_k
%
%   备注：需要配合SaveData.m一起使用，一共需要生成1000次估计数据    
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

%仿真步长
N = 1000;
LineWidth = 2;
h_con = sqrt(1.3);

q_con = 6;                   %系统噪声均值
r_con = 5;                   %测量噪声均值
Q_con = 5;                  %系统噪声方差矩阵
R_con = 0.25;                %测量噪声方差矩阵

q = q_con;                   %系统噪声均值
r = r_con;                   %测量噪声均值
Q = Q_con;                   %系统噪声方差矩阵
R = R_con;                   %测量噪声方差矩阵

%计算alpha阶次对应的GL定义系数 binomial coefficient 
bino_fir = zeros(1,N);       %微分阶次为0.7时GL定义下的系数
alpha = 0.7;
bino_fir(1,1) = 1;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);  
end

times = 100;
NumOrder = 1;
Q_point = zeros(1,times);
q_point = zeros(1,times);

for kk = 1:1:times
    
    NumOrder
    
    %状态测量初始化
    X_real = zeros(1,N);         %真实状态
    X_esti = zeros(1,N);         %状态最优估计值
    P_xesti = zeros(1,N);        %估计误差方差阵
    Z_meas = zeros(1,N);         %实际观测值

    %系统矩阵设置
    % A = [0,1; -0.1,-0.2];      %系统矩阵
    % B = [0; 1];                %
    % C = [0.1,0.3];             %
    I = eye(1,1);                %生成单位阵
    %I(3,3) = 0;

    %噪声
    W_noise = sqrt(Q)*randn(1,N) + q;  %系统噪声
    V_noise = sqrt(R)*randn(1,N) + r;  %测量噪声

    %初始值设置（初始矩阵不能为零）
    q_esti = zeros(1,N);         %系统噪声均值估计值
    r_esti = zeros(1,N);         %测量噪声均值估计值
    Q_esti = zeros(1,N);         %系统噪声方差矩阵估计值
    R_esti = zeros(1,N);         %测量噪声方差矩阵估计值

    P_pred_0 = 100;              %初始预测方差阵
    P_xesti(1,1) = P_pred_0;     %初始估计方差阵
    x_0  = 0;                    %初始状态     
    X_real(1,1) = x_0;           %真实状态初始值
    X_esti(1,1) = 0;             %状态估计初始值
    Z_meas(1,1) = V_noise(1,1);  %测量数据初始值

    %GL定义下短记忆原理的长度
    L = N+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % diff_X_real    表示k时刻状态的微分
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    diff_X_real = 0;

    %*************************************************************************%
    %     q_esti(k) = ( (k-1)*q_esti(k-1) + X_esti(1,k) - X_pre(k) + q_esti(k-1) )/k;
    %     Q_esti(k) = ( (k-1)*Q_esti(k-1) + (X_esti(1,k)-X_pre(k))^2 )/k;
    %     r_esti(k) = ( (k-1)*r_esti(k-1) + Z_meas(1,k) - Z_esti )/k;
    %     R_esti(k) = ( (k-1)*R_esti(k-1) + (Z_meas(1,k) - Z_esti)^2 )/k;
    %*************************************************************************%

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
            X_pre = diff_X_esti - rema + q;     %一步状态预测
            %预测误差协方差矩阵:P_pred

                %对误差矩阵进行cholsky分解
                S_chol = chol(P_xesti(1,k-1))';

                %计算余项
                rema_P = 0;
                if k>L+1
                    for i = 2:1:L+2
                        rema_P = rema_P + bino_fir(1,i)*P_xesti(1,k+1-i)*bino_fir(1,i)';
                    end
                else
                    for i = 2:1:k
                        rema_P = rema_P + bino_fir(1,i)*P_xesti(1,k+1-i)*bino_fir(1,i)';
                    end
                end

            %临时变量 temp_fun : 函数差值,函数为单变量
            temp_fun = 3*sin( 2*(X_esti(1,k-1)+h_con*S_chol) ) - (X_esti(1,k-1)+h_con*S_chol) - ...
                       3*sin( 2*(X_esti(1,k-1)-h_con*S_chol) ) + (X_esti(1,k-1)-h_con*S_chol);
            temp(k) = 1/(4*h_con^2) * temp_fun^2 + rema_P + ...
                      1/h_con*0.5*temp_fun*S_chol'*(-bino_fir(1,2))' + ...
                      1/h_con*(-bino_fir(1,2))*S_chol*0.5*temp_fun';
            P_xpred = temp(k) + Q_esti(k-1);
            %测量值估计  Z_esti ---- Z_k|k-1
            Z_esti = X_pre + r_esti(k-1);

            %测量预测误差协方差矩阵:P_zpred ---- P_z_k|k-1
            P_zpred = P_xpred + R;

            %计算卡尔曼增益:Kk(2*1)
            Kk = P_xpred/P_zpred;

            %r_esti(k) = ( (k-1)*r_esti(k-1) + Z_meas(1,k) - X_pre )/k;
            r_esti(k) = ( (k-1)*r_esti(k-1) + Z_meas(1,k) - X_pre )/k;

            %状态更新
            X_esti(1,k) = X_pre + Kk*( Z_meas(1,k) - Z_esti );

            %估计方差矩阵更新
            P_xesti(1,k) = P_zpred - Kk*P_zpred*Kk';

            Q_esti(k) = ( (k-1)*Q_esti(k-1) + ( X_esti(1,k) - X_pre )^2 - temp(k) + P_xesti(1,k) )/k;
            %q_esti(k) = ( (k-1)*q_esti(k-1) +   X_esti(1,k) - diff_X_esti + rema )/k;

    end

    Q_point(1,NumOrder) = Q_esti(k);
    r_point(1,NumOrder) = r_esti(k);

    NumOrder = NumOrder+1;

end

save estiqr Q_point r_point

% 前一百数据运行此段，之后直接运行SaveData.m
% Q_data = Q_point;
% r_data = r_point;
% save EstiData Q_data r_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%画图输出 均值方差估计散点图
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure;
% scatter(r_point(1,:),Q_point(1,:),'b','linewidth',LineWidth);
% legend('均值方差估计散点图','Location','best');
% hold on
% plot(r_con,Q_con,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],...
%     'MarkerSize',5, 'Marker','o','LineWidth',1, 'Color',[0 0 0]);
% hold off
% 
% figure;
% scatterhist(r_point,Q_point)
% hold on
% plot(r_con,Q_con,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],...
%     'MarkerSize',5, 'Marker','o','LineWidth',1, 'Color',[0 0 0]);
% hold off
% 
% figure;
% scatterhist(r_point,Q_point,'Kernel','on');
% 
