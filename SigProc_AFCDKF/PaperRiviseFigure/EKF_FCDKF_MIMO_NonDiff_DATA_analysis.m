%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   分数阶卡尔曼滤波器
%   论文：
%   备注：对获取的数据 (MIMO NonDiff FEKF VS FCDKF ) 进行比较
%         多输入多输出不完全可导情况下情况下FEKF与FCDKF获取的数据
%  
%   不可导状态分析测试
%           
%   备注: matlab中的数值精度比较高，很难出现状态 x_1 = 0 不可导的情况。因此
%        本次仿真中当状态 x_1 < 0.01 时，强行令 x_1 = 0 ，以检验算法的优越性
%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 50;
% k = 100;
% r = 1;
LineWidth = 1.5;


%获取状态数据文件，共三个状态三行，第一个状态会出现 0 的情况，即不可导的情况
 load EKF_FCDKF_MIMO_NonDiff_X_real            %X_real 保存的实际变量
 load EKF_FCDKF_MIMO_NonDiff_X_esti_FCDKF      %X_esti_FCDKF FCDKF估计的状态
 load EKF_FCDKF_MIMO_NonDiff_X_esti_FEKF       %X_esti_FEKF  FEKF估计的状态
 load EKF_FCDKF_MIMO_NonDiff_Outliers          %Outliers  异常值（不可导值）出现的位置
 
k = 1:1:N;

%状态 3 估计图用来验证FEKF的局限性以及FCDKF的有效性
figure;
%err_state = X_real - X_esti_A;
plot(k,X_real(3,:),'b',k,X_esti_FCDKF(3,:),'--r',k,X_esti_FEKF(3,:),'-.g', ...
    'linewidth',LineWidth);
% parameter_r = legend('$\hat{r}$','r','Location','best');
% set(parameter_r,'Interpreter','latex')
set(gcf,'Position',[200 200 400 300]); 
axis([0 35 -4 2]) %设置坐标轴在指定的区间
line([Outliers,Outliers],[-4,2],'linewidth',1,'color','k','LineStyle','--');
axis normal
ylabel('x_3','FontSize',7)
xlabel('iteration times','FontSize',8)
%设置坐标轴刻度字体名称，大小
set(gca,'FontName','Helvetica','FontSize',8)
%legend('real state','estimated state','Location','best');
legend('real value','FCDKF','FEKF','Location','best');
% 
% %所有状态估计图
% figure;
% %err_state = X_real - X_esti_A;
% plot(k,X_real(3,:),'b',k,X_esti_FCDKF(3,:),'--r',k,X_esti_FEKF(3,:),'-.g', ...
%     'linewidth',LineWidth);
% % parameter_r = legend('$\hat{r}$','r','Location','best');
% % set(parameter_r,'Interpreter','latex')
% set(gcf,'Position',[200 200 400 300]); 
% axis([0 25 -2 0.5]) %设置坐标轴在指定的区间
% line([Outliers,Outliers],[-2,0.5],'linewidth',1,'color','k','LineStyle','--');
% axis normal
% ylabel('x_3','FontSize',7)
% xlabel('iteration times','FontSize',8)
% %设置坐标轴刻度字体名称，大小
% set(gca,'FontName','Helvetica','FontSize',8)
%legend('real state','estimated state','Location','best');
legend('real value','FCDKF','FEKF','Location','best');