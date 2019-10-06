%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   �����׿������˲������渴��
%   ���ģ�     fractional order CDKF
%   Ŀ�ģ�FCDKF��FEKF�����ܱȽ�
%         ����ʵ��:    D^{0.7} x_k = 3*sin(2*x_{k-1}) -x_{k-1} + w_k
%                              y_k = x_k + v_k
%   �����
%
%   ��ע��
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

%���沽��
N = 200;

h_con = sqrt(3);

q_con = 1;                   %ϵͳ������ֵ
r_con = 1;                   %����������ֵ
Q_con = 0.81;                %ϵͳ�����������
R_con = 0.25;                %���������������

q = q_con;                   %ϵͳ������ֵ
r = r_con;                   %����������ֵ
Q = Q_con;                   %ϵͳ�����������
R = R_con;                   %���������������

%GL�����¶̼���ԭ��ĳ���
L = N+1;

%����alpha�״ζ�Ӧ��GL����ϵ�� binomial coefficient 
bino_fir = zeros(1,N);       %΢�ֽ״�Ϊ0.7ʱGL�����µ�ϵ��
alpha = 0.7;
bino_fir(1,1) = 1;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);  
end

%ϵͳ��������
% A = [0,1; -0.1,-0.2];      %ϵͳ����
% B = [0; 1];                %
% C = [0.1,0.3];             %
I = eye(1,1);                %���ɵ�λ��
%I(3,3) = 0;

%���ؿ����㷨
kk_N = 100;

err_state_FCDKF = zeros(kk_N,N);
err_state_FEKF  = zeros(kk_N,N);

for k_k = 1:1:kk_N

%״̬������ʼ��
X_real = zeros(1,N);         %��ʵ״̬
Z_meas = zeros(1,N);         %ʵ�ʹ۲�ֵ

%����
W_noise = sqrt(Q)*randn(1,N) + q;  %ϵͳ����
V_noise = sqrt(R)*randn(1,N) + r;  %��������

x_0  = 0;                    %��ʼ״̬     
X_real(1,1) = x_0;           %��ʵ״̬��ʼֵ
Z_meas(1,1) = V_noise(1,1);  %�������ݳ�ʼֵ
for k=2:1:N
    %����ʵ��״̬
    diff_X_real = 3*sin(2*X_real(1,k-1)) -X_real(1,k-1) + W_noise(1,k-1);
    rema = 0;
    for i = 2:1:k
        rema = rema + bino_fir(1,i)*X_real(1,k+1-i);
    end
    X_real(1,k) = diff_X_real - rema;
    %ʵ�ʹ۲�ֵ
    Z_meas(1,k) = X_real(1,k) + V_noise(1,k); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------�����׿������˲������ܲ���------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X_esti_A = zeros(1,N);       %״̬���Ź���ֵ
P_xesti = zeros(1,N);        %����������

P_pred_0 = 100;              %��ʼԤ�ⷽ����
P_xesti(1,1) = P_pred_0;     %��ʼ���Ʒ�����

for k=2:1:N
      %�������˲�
        %״̬Ԥ��:X_pre
        diff_X_esti = 3*sin(2*X_esti_A(1,k-1)) - X_esti_A(1,k-1);
            %��������
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
        X_pre = diff_X_esti - rema + q;     %һ��״̬Ԥ��
        %Ԥ�����Э�������:P_pred

            %�����������cholsky�ֽ�
            S_chol = chol(P_xesti(1,k-1))';

            %��������
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

        %��ʱ���� temp_fun : ������ֵ,����Ϊ������
        temp_fun = 3*sin( 2*(X_esti_A(1,k-1)+h_con*S_chol) ) - (X_esti_A(1,k-1)+h_con*S_chol) - ...
                   3*sin( 2*(X_esti_A(1,k-1)-h_con*S_chol) ) + (X_esti_A(1,k-1)-h_con*S_chol);
        temp = 1/(4*h_con^2) * temp_fun^2 + rema_P + ...
                  1/h_con*0.5*temp_fun*S_chol'*(-bino_fir(1,2))' + ...
                  1/h_con*(-bino_fir(1,2))*S_chol*0.5*temp_fun';
        P_xpred = temp + Q;
        %����ֵ����  Z_esti ---- Z_k|k-1
        Z_esti = X_pre + r;

        %����Ԥ�����Э�������:P_zpred ---- P_z_k|k-1
        P_zpred = P_xpred + R;

        %���㿨��������:Kk(2*1)
        Kk = P_xpred/P_zpred;

        %״̬����
        X_esti_A(1,k) = X_pre + Kk*( Z_meas(1,k) - Z_esti );

        %���Ʒ���������
        P_xesti(1,k) = P_zpred - Kk*P_zpred*Kk';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------��չ�����׿������˲������ܲ���----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X_esti_B = zeros(1,N);       %״̬���Ź���ֵ
P_xesti = zeros(1,N);        %����������

%��ʼֵ���ã���ʼ������Ϊ�㣩
P_pred_0 = 100;              %��ʼԤ�ⷽ����
P_xesti(1,1) = P_pred_0;     %��ʼ���Ʒ�����

for k=2:1:N
  %�������˲�
      %״̬Ԥ��:X_pre
        diff_X_esti = 3*sin(2*X_esti_B(1,k-1)) - X_esti_B(1,k-1);
            %��������
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
        X_pre = diff_X_esti - rema;     %һ��״̬Ԥ��
        %Ԥ�����Э�������:P_pred
            %��������
            rema_P = 0;
            if k>L+1
                for i = 3:1:L+2
                    rema_P = rema_P + bino_fir(1,i)*P_xesti(1,k+1-i)*bino_fir(1,i)';
                end
            else
                for i = 3:1:k
                    rema_P = rema_P + bino_fir(1,i)*P_xesti(1,k+1-i)*bino_fir(1,i)';
                end
            end
            F = 6*cos(2*X_esti_B(1,k-1)) - 1;
        P_xpred = (F-bino_fir(1,2))*P_xesti(1,k-1)*(F-bino_fir(1,2))'+ Q + rema_P;
        
        %���㿨��������:Kk(2*1)
        H = 1;
        Kk = P_xpred*H'/(H*P_xpred*H' + R);
        
        %״̬����
        X_esti_B(1,k) = X_pre + Kk*( Z_meas(1,k) - X_pre );
        
        %���Ʒ���������
        P_xesti(1,k) = (I-Kk*H)*P_xpred;
end

err_state_FCDKF(k_k,:) = (X_real - X_esti_A).^2; %FCDKF_ERROR
err_state_FEKF(k_k,:)  = (X_real - X_esti_B).^2; %FEKF_ERROR

end
%payphone/call me baby/wide awake/stars hip/we are young


% %״̬����ͼ
% k = 1:1:N;
% figure;
% plot(k,X_real(1,:),'b',k,X_esti_B(1,:),':r','linewidth',2);
% legend('ʵ��״̬','����״̬','Location','best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ͼ��� ��ֵ�������ɢ��ͼ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%������������ͼ
k = 1:1:N;

LineWidth = 1.5;

%ƽ�����ͼ
%square error
figure;
err1_state_FCDKF = (X_real - X_esti_A).^2;%FCDKF_ERROR
err1_state_FEKF = (X_real - X_esti_B).^2;%FEKF_ERROR
subplot(211)
plot(k,err1_state_FEKF(1,:),'b',k,err1_state_FCDKF(1,:),':r','linewidth',LineWidth);
%set(gcf,'Position',[200 200 400 300]); 
%axis([xmin xmax ymin ymax])������������ָ��������
axis normal
axis([ 0 N 0 6 ])
ylabel('SE','FontSize',8)
xlabel('time(sec)','FontSize',8)
%����������̶��������ƣ���С
set(gca,'FontName','Helvetica','FontSize',8)
legend('real state','estimated state','Location','best');
legend('FEKF','FCDKF','Location','best');

% RMSE ���������
% err_state_FCDKF = zeros(kk_N,N);
% err_state_FEKF  = zeros(kk_N,N);

% sum_FEKF  = 0;
% sum_FCDKF = 0;
RMSE_FEKF  = zeros(1,N); 
RMSE_FCDKF  = zeros(1,N);
for i = 1:1:N
    for j=1:1:kk_N
        RMSE_FEKF(1,i) =  RMSE_FEKF(1,i) + err_state_FEKF(j,i);
        RMSE_FCDKF(1,i) = RMSE_FCDKF(1,i) + err_state_FCDKF(j,i) ;
    end
end
%����RMSE
RMSE_FEKF = sqrt( RMSE_FEKF/kk_N );
RMSE_FCDKF = sqrt( RMSE_FCDKF/kk_N ) ;

k=1:1:N;
subplot(212)
plot(k,RMSE_FEKF(1,k),'b',k, RMSE_FCDKF(1,k),':r','linewidth',LineWidth);
%set(gcf,'Position',[200 200 400 300]); 
%axis([xmin xmax ymin ymax])������������ָ��������
axis normal
ylabel('RMSE','FontSize',8)
xlabel('time(sec)','FontSize',8)
%����������̶��������ƣ���С
set(gca,'FontName','Helvetica','FontSize',8)
%legend('FEKF','FCDKF','Location','best');
% figure_FontSize=8;
% set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
% set(findobj('FontSize',8.8),'FontSize',figure_FontSize);
set(gcf,'Position',[200 200 400 280]); 






