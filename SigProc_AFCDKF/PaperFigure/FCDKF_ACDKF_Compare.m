%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   分数阶卡尔曼滤波器仿真复现
%   论文：     fractional order CDKF
%   目的：噪声未知情况下，FCDKF与AFCDKF的性能比较
%         函数实验:    D^{0.7} x_k = 3*sin(2*x_{k-1}) -x_{k-1} + w_k
%                              y_k = x_k + v_k
%   结果：AFCDKF与FCDKF的性能比较
%
%   备注：
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

%仿真步长
N =10000;
LineWidth = 1.5;
h_con = sqrt(1.3);

q_con = 6;                   %系统噪声均值
r_con = 5;                   %测量噪声均值
Q_con = 10;                  %系统噪声方差矩阵
R_con = 0.25;                %测量噪声方差矩阵

q = q_con;                   %系统噪声均值
r = r_con;                   %测量噪声均值
Q = Q_con;                   %系统噪声方差矩阵
R = R_con;                   %测量噪声方差矩阵

%GL定义下短记忆原理的长度
L = N+1;

%计算alpha阶次对应的GL定义系数 binomial coefficient 
bino_fir = zeros(1,N);       %微分阶次为0.7时GL定义下的系数
alpha = 0.7;
bino_fir(1,1) = 1;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);  
end

%系统矩阵设置
% A = [0,1; -0.1,-0.2];      %系统矩阵
% B = [0; 1];                %
% C = [0.1,0.3];             %
I = eye(1,1);                %生成单位阵
%I(3,3) = 0;

%状态测量初始化
X_real = zeros(1,N);         %真实状态
Z_meas = zeros(1,N);         %实际观测值

%噪声
W_noise = sqrt(Q)*randn(1,N) + q;  %系统噪声
V_noise = sqrt(R)*randn(1,N) + r;  %测量噪声

x_0  = 0;                    %初始状态     
X_real(1,1) = x_0;           %真实状态初始值
Z_meas(1,1) = V_noise(1,1);  %测量数据初始值
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------自适应分数阶卡尔曼滤波器性能测试-----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X_esti_A = zeros(1,N);       %状态最优估计值
P_xesti = zeros(1,N);        %估计误差方差阵

%初始值设置（初始矩阵不能为零）
r_esti_A = zeros(1,N);       %测量噪声均值估计值
Q_esti_A = zeros(1,N);       %系统噪声方差矩阵估计值

P_pred_0 = 100;              %初始预测方差阵
P_xesti(1,1) = P_pred_0;     %初始估计方差阵

for k=2:1:N
      %卡尔曼滤波
        %状态预测:X_pre
        diff_X_esti = 3*sin(2*X_esti_A(1,k-1)) - X_esti_A(1,k-1);
            %计算余项
            rema = 0;
            if k>L
                for i = 2:1:L+1
                   rema = rema + bino_fir(1,i)*X_esti_A(1,k+1-i);
                end
            else
                for i = 2:1:k
                    rema = rema + bino_fir(1,i)*X_esti_A(1,k+1-i);
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
        temp_fun = 3*sin( 2*(X_esti_A(1,k-1)+h_con*S_chol) ) - (X_esti_A(1,k-1)+h_con*S_chol) - ...
                   3*sin( 2*(X_esti_A(1,k-1)-h_con*S_chol) ) + (X_esti_A(1,k-1)-h_con*S_chol);
        temp = 1/(4*h_con^2) * temp_fun^2 + rema_P + ...
                  1/h_con*0.5*temp_fun*S_chol'*(-bino_fir(1,2))' + ...
                  1/h_con*(-bino_fir(1,2))*S_chol*0.5*temp_fun';
        P_xpred = temp + Q_esti_A(k-1);
        %测量值估计  Z_esti ---- Z_k|k-1
        Z_esti = X_pre + r_esti_A(k-1);

        %测量预测误差协方差矩阵:P_zpred ---- P_z_k|k-1
        P_zpred = P_xpred + R;

        %计算卡尔曼增益:Kk(2*1)
        Kk = P_xpred/P_zpred;

        %r_esti(k) = ( (k-1)*r_esti(k-1) + Z_meas(1,k) - X_pre )/k;
        r_esti_A(k) = ( (k-1)*r_esti_A(k-1) + Z_meas(1,k) - X_pre )/k;

        %状态更新
        X_esti_A(1,k) = X_pre + Kk*( Z_meas(1,k) - Z_esti );

        %估计方差矩阵更新
        P_xesti(1,k) = P_zpred - Kk*P_zpred*Kk';

        Q_esti_A(k) = ( (k-1)*Q_esti_A(k-1) + ( X_esti_A(1,k) - X_pre )^2 - temp + P_xesti(1,k) )/k;
        %q_esti(k) = ( (k-1)*q_esti(k-1) +   X_esti(1,k) - diff_X_esti + rema )/k;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------普通分数阶卡尔曼滤波器性能测试-----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%未知系统噪声方差Q与测量噪声均值r
Q_guess = 8; %Q_real = 1
r_guess = 4; %r_real = 5
X_esti_B = zeros(1,N);       %状态最优估计值
P_xesti = zeros(1,N);        %估计误差方差阵

%初始值设置（初始矩阵不能为零）
P_pred_0 = 100;              %初始预测方差阵
P_xesti(1,1) = P_pred_0;     %初始估计方差阵

for k=2:1:N
      %卡尔曼滤波
        %状态预测:X_pre
        diff_X_esti = 3*sin(2*X_esti_B(1,k-1)) - X_esti_B(1,k-1);
            %计算余项
            rema = 0;
            if k>L
                for i = 2:1:L+1
                   rema = rema + bino_fir(1,i)*X_esti_B(1,k+1-i);
                end
            else
                for i = 2:1:k
                    rema = rema + bino_fir(1,i)*X_esti_B(1,k+1-i);
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
        temp_fun = 3*sin( 2*(X_esti_B(1,k-1)+h_con*S_chol) ) - (X_esti_B(1,k-1)+h_con*S_chol) - ...
                   3*sin( 2*(X_esti_B(1,k-1)-h_con*S_chol) ) + (X_esti_B(1,k-1)-h_con*S_chol);
        temp = 1/(4*h_con^2) * temp_fun^2 + rema_P + ...
                  1/h_con*0.5*temp_fun*S_chol'*(-bino_fir(1,2))' + ...
                  1/h_con*(-bino_fir(1,2))*S_chol*0.5*temp_fun';
        P_xpred = temp + Q_guess;
        %测量值估计  Z_esti ---- Z_k|k-1
        Z_esti = X_pre + r_guess;

        %测量预测误差协方差矩阵:P_zpred ---- P_z_k|k-1
        P_zpred = P_xpred + R;

        %计算卡尔曼增益:Kk(2*1)
        Kk = P_xpred/P_zpred;

        %状态更新
        X_esti_B(1,k) = X_pre + Kk*( Z_meas(1,k) - Z_esti );

        %估计方差矩阵更新
        P_xesti(1,k) = P_zpred - Kk*P_zpred*Kk';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%画图输出 均值方差估计散点图
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%输入与测量输出图
% figure;
% plot(k,Z_meas,':r',k,X_real(1,:),'b','linewidth',LineWidth);
% legend('测量输出','Location','best');

%状态估计图
% figure;
% plot(k,X_real(1,:),'b',k,X_esti_A(1,:),':y',k,X_esti_B(1,:),':r','linewidth',1);
% legend('实际状态','估计状态','Location','best');

%误差图
%误差计算
% err_state = zeros(2,N);

k = 1:1:N;
figure;

plot(k,Q_esti_A(1,:),'b',k,r_esti_A(1,:),'r','linewidth',LineWidth);
axis([0 1000 0 25])%设置坐标轴在指定的区间
set(gcf,'Position',[200 200 400 300]); 
%axis([xmin xmax ymin ymax])设置坐标轴在指定的区间
axis normal
ylabel('estimated value','FontSize',8)
xlabel('iteration times','FontSize',8)
%设置坐标轴刻度字体名称，大小
set(gca,'FontName','Helvetica','FontSize',8)
line([0 1000],[r r],'Color','r','LineStyle','--')
line([0 1000],[Q Q],'Color','b','LineStyle','--')
legend('$\hat{r}$','$\hat{Q}$','Location','best');

k = 1:1:N;
figure;
err_state_A = (X_real - X_esti_A).^2;%AFCDKF_ERROR
err_state_B = (X_real - X_esti_B).^2;%FCDKF_ERROR
plot(k,err_state_B(1,:),'b',k,err_state_A(1,:),'g','linewidth',LineWidth);
axis([0 N 0 12])%设置坐标轴在指定的区间
set(gcf,'Position',[200 200 400 300]); 
%axis([xmin xmax ymin ymax])设置坐标轴在指定的区间
axis normal
ylabel('SE','FontSize',8)
xlabel('iteration times','FontSize',8)
%设置坐标轴刻度字体名称，大小
set(gca,'FontName','Helvetica','FontSize',8)
legend('FCDKF','AFCDKF','Location','best');

