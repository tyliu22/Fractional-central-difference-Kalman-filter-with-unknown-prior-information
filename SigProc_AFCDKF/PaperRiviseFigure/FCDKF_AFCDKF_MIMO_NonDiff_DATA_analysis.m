%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   �����׿������˲���
%   ���ģ�
%   ��ע���Ի�ȡ������ (MIMO NonDiff AFCDKF VS FCDKF ) ���бȽ�
%         ��������������ȫ�ɵ�����������FCDKF��AFCDKF��ȡ������
%         �����м���������ֵδ֪
%  
%           
%  ����
%
%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 50;
% k = 100;
% r = 1;
LineWidth = 1.5;
k = 1:1:N;

%%
% save FCDKF_AFCDKF_MIMO_NonDiff_X_real        X_real
% save FCDKF_AFCDKF_MIMO_NonDiff_X_esti_AFCDKF X_esti_AFCDKF
% save FCDKF_AFCDKF_MIMO_NonDiff_X_esti_FCDKF  X_esti_FCDKF
% save FCDKF_AFCDKF_MIMO_NonDiff_r_esti        r_esti
% save FCDKF_AFCDKF_MIMO_NonDiff_r_real        r_real
%��ȡ״̬�����ļ���������״̬���У����ڷ��������㷨������
 load FCDKF_AFCDKF_MIMO_NonDiff_X_real         % X_real         �����ʵ�ʱ���
 load FCDKF_AFCDKF_MIMO_NonDiff_X_esti_AFCDKF  % X_esti_AFCDKF  AFCDKF���Ƶ�״̬
 load FCDKF_AFCDKF_MIMO_NonDiff_X_esti_FCDKF   % X_esti_FCDKF   FCDKF���Ƶ�״̬
 load FCDKF_AFCDKF_MIMO_NonDiff_r_esti         % r_esti         AFCDKF���ƵĲ���������ֵ
 load FCDKF_AFCDKF_MIMO_NonDiff_r_real         % r_real         ʵ�ʲ��������ľ�ֵ

 % ���� SE ƽ�����
 err_state_AFCDKF1 = ( X_real(1,:) - X_esti_AFCDKF(1,:) ).^2; %
 err_state_FCDKF1  = ( X_real(1,:) - X_esti_FCDKF(1,:)  ).^2; %
 err_state_AFCDKF2 = ( X_real(2,:) - X_esti_AFCDKF(2,:) ).^2; %
 err_state_FCDKF2  = ( X_real(2,:) - X_esti_FCDKF(2,:)  ).^2; %
 err_state_AFCDKF3 = ( X_real(3,:) - X_esti_AFCDKF(3,:) ).^2; %
 err_state_FCDKF3  = ( X_real(3,:) - X_esti_FCDKF(3,:)  ).^2; %
 
 %% RMSE ������
 RMSE_AFCDKF1 = zeros(1,N);
 RMSE_FCDKF1  = zeros(1,N);
 RMSE_AFCDKF2 = zeros(1,N);
 RMSE_FCDKF2  = zeros(1,N);
 RMSE_AFCDKF3 = zeros(1,N);
 RMSE_FCDKF3  = zeros(1,N);
 
     RMSE_AFCDKF1(1,1) = err_state_AFCDKF1(1,1);
     RMSE_FCDKF1(1,1)  = err_state_FCDKF1(1,1) ;
     RMSE_AFCDKF2(1,1) = err_state_AFCDKF2(1,1);
     RMSE_FCDKF2(1,1)  = err_state_FCDKF2(1,1) ;
     RMSE_AFCDKF3(1,1) = err_state_AFCDKF3(1,1);
     RMSE_FCDKF3(1,1)  = err_state_FCDKF3(1,1) ;
     
 for i = 2:1:N
     RMSE_AFCDKF1(1,i) = RMSE_AFCDKF1(1,i-1) + err_state_AFCDKF1(1,i);
     RMSE_FCDKF1(1,i)  = RMSE_FCDKF1(1,i-1)  + err_state_FCDKF1(1,i) ;
     RMSE_AFCDKF2(1,i) = RMSE_AFCDKF2(1,i-1) + err_state_AFCDKF2(1,i);
     RMSE_FCDKF2(1,i)  = RMSE_FCDKF2(1,i-1)  + err_state_FCDKF2(1,i) ;
     RMSE_AFCDKF3(1,i) = RMSE_AFCDKF3(1,i-1) + err_state_AFCDKF3(1,i);
     RMSE_FCDKF3(1,i)  = RMSE_FCDKF3(1,i-1)  + err_state_FCDKF3(1,i) ;
 end
%����RMSE
 for i = 1:1:N
     RMSE_AFCDKF1(1,i) = sqrt( RMSE_AFCDKF1(1,i) / i );
     RMSE_FCDKF1(1,i)  = sqrt( RMSE_FCDKF1(1,i)  / i );
     RMSE_AFCDKF2(1,i) = sqrt( RMSE_AFCDKF2(1,i) / i );
     RMSE_FCDKF2(1,i)  = sqrt( RMSE_FCDKF2(1,i)  / i );
     RMSE_AFCDKF3(1,i) = sqrt( RMSE_AFCDKF3(1,i) / i );
     RMSE_FCDKF3(1,i)  = sqrt( RMSE_FCDKF3(1,i)  / i );
 end
 
 figure;
set(gcf,'Position',[200 200 400 300]); 

subplot(311)
plot(k,RMSE_AFCDKF1(1,:),'b',k,RMSE_FCDKF1(1,:),'--r','linewidth',LineWidth);
%set(gcf,'Position',[200 200 400 300]); 
%axis([0 N 0 30]) %axis([xmin xmax ymin ymax])������������ָ��������
axis normal
ylabel('RMSE(x_1)','FontSize',7)
xlabel('iteration times','FontSize',8)
%����������̶��������ƣ���С
%legend('real state','estimated state','Location','best');
%legend('FCDKF','FEKF','Location','best');

subplot(312)
plot(k,RMSE_AFCDKF2(1,:),'b',k,RMSE_FCDKF2(1,:),'--r','linewidth',LineWidth);
%set(gcf,'Position',[200 200 400 300]); 
%axis([0 N 0 5]) 
%axis([xmin xmax ymin ymax])������������ָ��������
axis normal
ylabel('RMSE(x_2)','FontSize',7)
xlabel('iteration times','FontSize',8)
%����������̶��������ƣ���С
set(gca,'FontName','Helvetica','FontSize',8)
% legend('real state','estimated state','Location','best');
% legend('FCDKF','FEKF','Location','best');

subplot(313)
plot(k,RMSE_AFCDKF3(1,:),'b',k,RMSE_FCDKF3(1,:),'--r','linewidth',LineWidth);
%set(gcf,'Position',[200 200 400 300]); 
%axis([0 N 0 0.02]) 
%axis([xmin xmax ymin ymax])������������ָ��������
axis normal
ylabel('RMSE(x_3)','FontSize',7)
xlabel('iteration times','FontSize',8)
%����������̶��������ƣ���С
set(gca,'FontName','Helvetica','FontSize',8)
% legend('real state','estimated state','Location','best');
 legend('AFCDKF','FCDKF','Location','best');

 
%%  ����������ֵ r ��������ͼ
 figure;
 plot(k,r_esti(1,k),'b','linewidth',LineWidth);
 line([0,N],[r_real,r_real],'linewidth',1.5,'color','r','LineStyle','--');
 parameter_r = legend('$\hat{r}$','r','Location','best');
 set(parameter_r,'Interpreter','latex')
 set(gcf,'Position',[200 200 400 300]); 
 axis([0 N 0.5 1.1]) %axis([xmin xmax ymin ymax])������������ָ��������
 axis normal
 ylabel('estimated value','FontSize',7)
 xlabel('iteration times','FontSize',8)
 %����������̶��������ƣ���С
 set(gca,'FontName','Helvetica','FontSize',8)
 %legend('real state','estimated state','Location','best');
 %legend('FCDKF','FEKF','Location','best');
 %err_state = X_real - X_esti_A;

%% SE ƽ�����ͼ
 %square error
 figure;
 set(gcf,'Position',[200 200 400 300]); 

 subplot(311)
 plot(k,err_state_AFCDKF1(1,:),'b',k,err_state_FCDKF1(1,:),'--r','linewidth',LineWidth);
 %set(gcf,'Position',[200 200 400 300]); 
 %axis([0 N 0 30]) %axis([xmin xmax ymin ymax])������������ָ��������
 axis normal
 ylabel('SE(x_1)','FontSize',7)
 xlabel('iteration times','FontSize',8)
 %����������̶��������ƣ���С
 %legend('real state','estimated state','Location','best');
 %legend('FCDKF','FEKF','Location','best');

 subplot(312)
 plot(k,err_state_AFCDKF2(1,:),'b',k,err_state_FCDKF2(1,:),'--r','linewidth',LineWidth);
 %set(gcf,'Position',[200 200 400 300]); 
 %axis([0 N 0 5]) 
 %axis([xmin xmax ymin ymax])������������ָ��������
 axis normal
 ylabel('SE(x_2)','FontSize',7)
 xlabel('iteration times','FontSize',8)
 %����������̶��������ƣ���С
 set(gca,'FontName','Helvetica','FontSize',8)
 % legend('real state','estimated state','Location','best');
 % legend('FCDKF','FEKF','Location','best');

 subplot(313)
 plot(k,err_state_AFCDKF3(1,:),'b',k,err_state_FCDKF3(1,:),'--r','linewidth',LineWidth);
 %set(gcf,'Position',[200 200 400 300]); 
 %axis([0 N 0 0.02]) 
 %axis([xmin xmax ymin ymax])������������ָ��������
 axis normal
 ylabel('SE(x_3)','FontSize',7)
 xlabel('iteration times','FontSize',8)
 %����������̶��������ƣ���С
 set(gca,'FontName','Helvetica','FontSize',8)
 % legend('real state','estimated state','Location','best');
 legend('AFCDKF','FCDKF','Location','best');


