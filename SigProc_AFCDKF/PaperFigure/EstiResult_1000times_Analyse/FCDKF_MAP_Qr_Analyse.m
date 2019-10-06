%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   自适应分数阶中心差分卡尔曼滤波器估计性能分析
%   论文：     fractional order CDKF
%   目的：将生成的1000次估计数据进行数据分析
%         函数实验:    D^{0.7} x_k = 3*sin(2*x_{k-1}) -x_{k-1} + w_k
%                              y_k = x_k + v_k
%
%   备注：    
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear

load EstiData.mat

r_con = 5;                  %测量噪声均值
Q_con = 5;                  %系统噪声方差矩阵

figure;
scatter(r_data(1,:),Q_data(1,:),'b','linewidth',1.5);
legend('均值方差估计散点图','Location','best');
hold on
plot(r_con,Q_con,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],...
    'MarkerSize',5, 'Marker','o','LineWidth',1, 'Color',[0 0 0]);
hold off
set(gcf,'Position',[200 200 400 300]); 
%axis([xmin xmax ymin ymax])设置坐标轴在指定的区间
axis normal
ylabel('covariance','FontSize',8)
xlabel('mean','FontSize',8)
%设置坐标轴刻度字体名称，大小
set(gca,'FontName','Helvetica','FontSize',8)
legend('estimated value','Location','best');



figure;
fig_sh = scatterhist(r_data,Q_data,'MarkerSize',3,'Marker','o');%'NBins',[30 30]
hold on
plot(r_con,Q_con,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],...
    'MarkerSize',5, 'Marker','o','LineWidth',1, 'Color',[0 0 0]);
hold off
set(gcf,'Position',[200 200 400 300]); 
%axis([xmin xmax ymin ymax])设置坐标轴在指定的区间
axis normal
ylabel('$$\hat{Q}$$','FontSize',8)
xlabel('$$\hat{r}$$','FontSize',8)
%设置坐标轴刻度字体名称，大小
set(gca,'FontName','Helvetica','FontSize',8)
legend('estimated value','Location','best');


% 
% set(gcf,'Position',[100 100 260 220]);%设置图像大小，可直接导为eps格式
% set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
% set(findobj('FontSize',10),'FontSize',figure_FontSize);%
% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);%设置线宽为2
% % figure;
% % scatterhist(r_data,Q_data,'Kernel','on');
% % hold on
% % plot(r_con,Q_con,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],...
% %     'MarkerSize',5, 'Marker','o','LineWidth',1, 'Color',[0 0 0]);
% % hold off