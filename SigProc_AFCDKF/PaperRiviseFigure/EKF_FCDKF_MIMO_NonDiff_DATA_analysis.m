%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   �����׿������˲���
%   ���ģ�
%   ��ע���Ի�ȡ������ (MIMO NonDiff FEKF VS FCDKF ) ���бȽ�
%         ��������������ȫ�ɵ�����������FEKF��FCDKF��ȡ������
%  
%   ���ɵ�״̬��������
%           
%   ��ע: matlab�е���ֵ���ȱȽϸߣ����ѳ���״̬ x_1 = 0 ���ɵ�����������
%        ���η����е�״̬ x_1 < 0.01 ʱ��ǿ���� x_1 = 0 ���Լ����㷨����Խ��
%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 50;
% k = 100;
% r = 1;
LineWidth = 1.5;


%��ȡ״̬�����ļ���������״̬���У���һ��״̬����� 0 ������������ɵ������
 load EKF_FCDKF_MIMO_NonDiff_X_real            %X_real �����ʵ�ʱ���
 load EKF_FCDKF_MIMO_NonDiff_X_esti_FCDKF      %X_esti_FCDKF FCDKF���Ƶ�״̬
 load EKF_FCDKF_MIMO_NonDiff_X_esti_FEKF       %X_esti_FEKF  FEKF���Ƶ�״̬
 load EKF_FCDKF_MIMO_NonDiff_Outliers          %Outliers  �쳣ֵ�����ɵ�ֵ�����ֵ�λ��
 
k = 1:1:N;

%״̬ 3 ����ͼ������֤FEKF�ľ������Լ�FCDKF����Ч��
figure;
%err_state = X_real - X_esti_A;
plot(k,X_real(3,:),'b',k,X_esti_FCDKF(3,:),'--r',k,X_esti_FEKF(3,:),'-.g', ...
    'linewidth',LineWidth);
% parameter_r = legend('$\hat{r}$','r','Location','best');
% set(parameter_r,'Interpreter','latex')
set(gcf,'Position',[200 200 400 300]); 
axis([0 35 -4 2]) %������������ָ��������
line([Outliers,Outliers],[-4,2],'linewidth',1,'color','k','LineStyle','--');
axis normal
ylabel('x_3','FontSize',7)
xlabel('iteration times','FontSize',8)
%����������̶��������ƣ���С
set(gca,'FontName','Helvetica','FontSize',8)
%legend('real state','estimated state','Location','best');
legend('real value','FCDKF','FEKF','Location','best');
% 
% %����״̬����ͼ
% figure;
% %err_state = X_real - X_esti_A;
% plot(k,X_real(3,:),'b',k,X_esti_FCDKF(3,:),'--r',k,X_esti_FEKF(3,:),'-.g', ...
%     'linewidth',LineWidth);
% % parameter_r = legend('$\hat{r}$','r','Location','best');
% % set(parameter_r,'Interpreter','latex')
% set(gcf,'Position',[200 200 400 300]); 
% axis([0 25 -2 0.5]) %������������ָ��������
% line([Outliers,Outliers],[-2,0.5],'linewidth',1,'color','k','LineStyle','--');
% axis normal
% ylabel('x_3','FontSize',7)
% xlabel('iteration times','FontSize',8)
% %����������̶��������ƣ���С
% set(gca,'FontName','Helvetica','FontSize',8)
%legend('real state','estimated state','Location','best');
legend('real value','FCDKF','FEKF','Location','best');