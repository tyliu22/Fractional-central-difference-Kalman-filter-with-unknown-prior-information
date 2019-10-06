%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   �����׿������˲������渴��
%   ���ģ�SP_AFCDKF
%   ��ע��FEKF VS AFCDKF ���ݷ��� 
%         ϵͳ�����ɵ����
%         MIMO����¼���δ֪����������ֵFEKF��AFCDKF�ıȽ�
%         �����ź�Ϊϵͳ��ʶ�г��õĶ������źţ��ڸ���������£���ֵ�ı�ʶ
%         Ч���Ϻ�
%   ʵ��״̬����
%           
%   ��ע��
%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
N = 100;
k = 100;
r = 1;
LineWidth = 2;

load FEKF_AFCDKF_MIMO_COMPARISON_DATA
    %X_real    ��ʾϵͳ��ʵ��״̬
    %X_esti_A  ��ʾAFCDKF�㷨���Ƶ�ϵͳ״̬
    %X_esti_B  ��ʾFEKF�㷨���Ƶ�ϵͳ״̬
    %r_esti_A  ��ʾAFCDKF�㷨���ƵĲ���������ֵ

k = 1:1:N;

%��������ͼ
figure;
err_state = X_real - X_esti_A;
plot(k,r_esti_A(1,k),'b','linewidth',LineWidth);
line([0,N],[r,r],'linewidth',1,'color','b','LineStyle','--');
parameter_r = legend('$\hat{r}$','r','Location','best');
set(parameter_r,'Interpreter','latex')
set(gcf,'Position',[200 200 400 300]); 
%axis([xmin xmax ymin ymax])������������ָ��������
axis normal
ylabel('estimated value','FontSize',7)
xlabel('iteration times','FontSize',8)
%����������̶��������ƣ���С
set(gca,'FontName','Helvetica','FontSize',8)
%legend('real state','estimated state','Location','best');
%legend('FCDKF','FEKF','Location','best');




%ƽ�����ͼ
%square error
figure;
set(gcf,'Position',[200 200 400 300]); 
err_state_A1 = ( X_real(1,:) - X_esti_A(1,:) ).^2; %
err_state_B1 = ( X_real(1,:) - X_esti_B(1,:) ).^2; %
err_state_A2 = ( X_real(2,:) - X_esti_A(2,:) ).^2; %
err_state_B2 = ( X_real(2,:) - X_esti_B(2,:) ).^2; %
err_state_A3 = ( X_real(3,:) - X_esti_A(3,:) ).^2; %
err_state_B3 = ( X_real(3,:) - X_esti_B(3,:) ).^2; %

subplot(311)
plot(k,err_state_A1(1,:),'b',k,err_state_B1(1,:),'--r','linewidth',LineWidth);
%set(gcf,'Position',[200 200 400 300]); 
axis([0 100 0 30]) %axis([xmin xmax ymin ymax])������������ָ��������
axis normal
ylabel('SE(x_1)','FontSize',7)
xlabel('iteration times','FontSize',8)
%����������̶��������ƣ���С
%legend('real state','estimated state','Location','best');
%legend('FCDKF','FEKF','Location','best');

subplot(312)
plot(k,err_state_A2(1,:),'b',k,err_state_B2(1,:),'--r','linewidth',LineWidth);
%set(gcf,'Position',[200 200 400 300]); 
axis([0 100 0 5]) 
%axis([xmin xmax ymin ymax])������������ָ��������
axis normal
ylabel('SE(x_2)','FontSize',7)
xlabel('iteration times','FontSize',8)
%����������̶��������ƣ���С
set(gca,'FontName','Helvetica','FontSize',8)
% legend('real state','estimated state','Location','best');
% legend('FCDKF','FEKF','Location','best');

subplot(313)
plot(k,err_state_A3(1,:),'b',k,err_state_B3(1,:),'--r','linewidth',LineWidth);
%set(gcf,'Position',[200 200 400 300]); 
axis([0 100 0 0.02]) 
%axis([xmin xmax ymin ymax])������������ָ��������
axis normal
ylabel('SE(x_3)','FontSize',7)
xlabel('iteration times','FontSize',8)
%����������̶��������ƣ���С
set(gca,'FontName','Helvetica','FontSize',8)
% legend('real state','estimated state','Location','best');
 legend('AFCDKF','FEKF','Location','best');



%״̬����ͼ
% k = 1:1:N;
% figure;
% plot(k,X_real(3,:),'b',k,X_esti_A(3,:),'--r',k,X_esti_B(3,:),':g','linewidth',LineWidth);
% legend('ʵ��״̬3','FCDKF����״̬3','FEKF����״̬3','Location','best');
% 
% figure;
% plot(k,X_real(2,:),'b',k,X_esti_A(2,:),'--r',k,X_esti_B(2,:),':g','linewidth',LineWidth);
% legend('ʵ��״̬2','FCDKF����״̬2','FEKF����״̬2','Location','best');
% 
% hold on
% plot(k,X_real(1,:),'b',k,X_esti_A(1,:),'--r',k,X_esti_B(1,:),':g','linewidth',LineWidth);
% legend('ʵ��״̬1','FCDKF����״̬1','FEKF����״̬1','Location','best');





