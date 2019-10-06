%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   �����׿������˲������渴��
%   ���ģ�     fractional order CDKF
%   Ŀ�ģ�����δ֪����£�FCDKF��AFCDKF�����ܱȽ�
%         ����ʵ��:    D^{0.7} x_k = 3*sin(2*x_{k-1}) -x_{k-1} + w_k
%                              y_k = x_k + v_k
%   �����AFCDKF��FCDKF�����ܱȽ�
%
%   ��ע��
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

%���沽��
N =10000;
LineWidth = 1.5;
h_con = sqrt(1.3);

q_con = 6;                   %ϵͳ������ֵ
r_con = 5;                   %����������ֵ
Q_con = 10;                  %ϵͳ�����������
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
%---------------------����Ӧ�����׿������˲������ܲ���-----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X_esti_A = zeros(1,N);       %״̬���Ź���ֵ
P_xesti = zeros(1,N);        %����������

%��ʼֵ���ã���ʼ������Ϊ�㣩
r_esti_A = zeros(1,N);       %����������ֵ����ֵ
Q_esti_A = zeros(1,N);       %ϵͳ��������������ֵ

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
        P_xpred = temp + Q_esti_A(k-1);
        %����ֵ����  Z_esti ---- Z_k|k-1
        Z_esti = X_pre + r_esti_A(k-1);

        %����Ԥ�����Э�������:P_zpred ---- P_z_k|k-1
        P_zpred = P_xpred + R;

        %���㿨��������:Kk(2*1)
        Kk = P_xpred/P_zpred;

        %r_esti(k) = ( (k-1)*r_esti(k-1) + Z_meas(1,k) - X_pre )/k;
        r_esti_A(k) = ( (k-1)*r_esti_A(k-1) + Z_meas(1,k) - X_pre )/k;

        %״̬����
        X_esti_A(1,k) = X_pre + Kk*( Z_meas(1,k) - Z_esti );

        %���Ʒ���������
        P_xesti(1,k) = P_zpred - Kk*P_zpred*Kk';

        Q_esti_A(k) = ( (k-1)*Q_esti_A(k-1) + ( X_esti_A(1,k) - X_pre )^2 - temp + P_xesti(1,k) )/k;
        %q_esti(k) = ( (k-1)*q_esti(k-1) +   X_esti(1,k) - diff_X_esti + rema )/k;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------��ͨ�����׿������˲������ܲ���-----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%δ֪ϵͳ��������Q�����������ֵr
Q_guess = 8; %Q_real = 1
r_guess = 4; %r_real = 5
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
        temp_fun = 3*sin( 2*(X_esti_B(1,k-1)+h_con*S_chol) ) - (X_esti_B(1,k-1)+h_con*S_chol) - ...
                   3*sin( 2*(X_esti_B(1,k-1)-h_con*S_chol) ) + (X_esti_B(1,k-1)-h_con*S_chol);
        temp = 1/(4*h_con^2) * temp_fun^2 + rema_P + ...
                  1/h_con*0.5*temp_fun*S_chol'*(-bino_fir(1,2))' + ...
                  1/h_con*(-bino_fir(1,2))*S_chol*0.5*temp_fun';
        P_xpred = temp + Q_guess;
        %����ֵ����  Z_esti ---- Z_k|k-1
        Z_esti = X_pre + r_guess;

        %����Ԥ�����Э�������:P_zpred ---- P_z_k|k-1
        P_zpred = P_xpred + R;

        %���㿨��������:Kk(2*1)
        Kk = P_xpred/P_zpred;

        %״̬����
        X_esti_B(1,k) = X_pre + Kk*( Z_meas(1,k) - Z_esti );

        %���Ʒ���������
        P_xesti(1,k) = P_zpred - Kk*P_zpred*Kk';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ͼ��� ��ֵ�������ɢ��ͼ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%������������ͼ
% figure;
% plot(k,Z_meas,':r',k,X_real(1,:),'b','linewidth',LineWidth);
% legend('�������','Location','best');

%״̬����ͼ
% figure;
% plot(k,X_real(1,:),'b',k,X_esti_A(1,:),':y',k,X_esti_B(1,:),':r','linewidth',1);
% legend('ʵ��״̬','����״̬','Location','best');

%���ͼ
%������
% err_state = zeros(2,N);

k = 1:1:N;
figure;

plot(k,Q_esti_A(1,:),'b',k,r_esti_A(1,:),'r','linewidth',LineWidth);
axis([0 1000 0 25])%������������ָ��������
set(gcf,'Position',[200 200 400 300]); 
%axis([xmin xmax ymin ymax])������������ָ��������
axis normal
ylabel('estimated value','FontSize',8)
xlabel('iteration times','FontSize',8)
%����������̶��������ƣ���С
set(gca,'FontName','Helvetica','FontSize',8)
line([0 1000],[r r],'Color','r','LineStyle','--')
line([0 1000],[Q Q],'Color','b','LineStyle','--')
legend('$\hat{r}$','$\hat{Q}$','Location','best');

k = 1:1:N;
figure;
err_state_A = (X_real - X_esti_A).^2;%AFCDKF_ERROR
err_state_B = (X_real - X_esti_B).^2;%FCDKF_ERROR
plot(k,err_state_B(1,:),'b',k,err_state_A(1,:),'g','linewidth',LineWidth);
axis([0 N 0 12])%������������ָ��������
set(gcf,'Position',[200 200 400 300]); 
%axis([xmin xmax ymin ymax])������������ָ��������
axis normal
ylabel('SE','FontSize',8)
xlabel('iteration times','FontSize',8)
%����������̶��������ƣ���С
set(gca,'FontName','Helvetica','FontSize',8)
legend('FCDKF','AFCDKF','Location','best');

