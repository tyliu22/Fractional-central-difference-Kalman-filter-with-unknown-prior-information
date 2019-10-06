%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   �����׿������˲������渴��
%   ���ģ�SP_AFCDKF
%   ��ע��FEKF VS AFCDKF 
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

%���沽��
N = 100;
LineWidth = 1.5;
h_con = sqrt(3);

%״̬������ʼ��
X_real = zeros(3,N);            %��ʵ״̬
Z_meas = zeros(1,N);            %ʵ�ʹ۲�ֵ

%ϵͳ��������
% A = [0,1; -0.1,-0.2];           %ϵͳ����
% B = [0; 1];                     %
% C = [0.1,0.3];                  %
I = eye(3,3);                     %���ɵ�λ��
%I(3,3) = 0;

%����
q = [0,0,0]';                      %ϵͳ������ֵ
r = 1;                             %����������ֵ
Q = [0.3,0,0; 0,0.3,0; 0,0,0.001]; %ϵͳ�����������
R = 0.3;                           %���������������
W_noise = Q*randn(3,N) + q;        %ϵͳ����
V_noise = R*randn(1,N) + r;        %��������

%��ʼֵ����
P_pred_0 = [100,0,0; 0,100,0;0,0,100];%��ʼԤ�ⷽ����
x_0  = [0; 0; 0.2];                   %��ʼ״̬
X_real(:,1) = x_0;                    %��ʵ״̬��ʼֵ
Z_meas(:,1) = V_noise(:,1);           %�������ݳ�ʼֵ

%����alpha�״ζ�Ӧ��GL����ϵ�� binomial coefficient 
bino_fir = zeros(1,N);          %΢�ֽ״�Ϊ0.7ʱGL�����µ�ϵ��
alpha = 0.7;
bino_fir(1,1) = 0.7;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);  
end

bino_sec = zeros(1,N);          %΢�ֽ״�Ϊ1.2ʱGL�����µ�ϵ��
alpha = 1.2;
bino_sec(1,1) = 1.2;
for i = 2:1:N
    bino_sec(1,i) = (1-(alpha+1)/(i-1))*bino_sec(1,i-1);  
end

bino_thi = zeros(1,N);          %΢�ֽ״�Ϊ1ʱGL�����µ�ϵ��
alpha = 1;
bino_thi(1,1) = 1;
for i = 2:1:N
    bino_thi(1,i) = (1-(alpha+1)/(i-1))*bino_thi(1,i-1);  
end

%����GL΢�ֶ�����ϵ������
gamma = cell(1,N);
temp_matrx = zeros(3,3);
for i = 1:1:N 
    temp_matrx(1,1) = bino_fir(1,i);
    temp_matrx(2,2) = bino_sec(1,i);
    temp_matrx(3,3) = bino_thi(1,i);
    gamma{1,i} = temp_matrx;
end

%ϵͳ��������:ǰһ������Ϊ1����һ������Ϊ

U_input = idinput(100,'sine',[0,1],[-1,1])';
% 'rgs'��˹����źţ�
% 'rbs'��Ĭ�ϣ���ֵ����ź� ��
% 'prbs'��ֵα����źţ�M���У�
% 'sine'�����źź�

% U_input = ones(1,N);            
% for i = N/2+1:1:N
%     U_input(:,i) = -2;
% end

%GL�����¶̼���ԭ��ĳ���
L = 20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diff_X_real    ��ʾkʱ��״̬��΢��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff_X_real = [0;0;0];

for k=2:1:N
    %����ʵ��״̬
    diff_X_real = [X_real(2,k-1); ...
                  -0.1*X_real(1,k-1)-X_real(2,k-1)*X_real(3,k-1)+ ...
                  U_input(1,k-1); 0] + W_noise(:,k-1);
    rema = [0;0;0];
    for i = 2:1:k
        rema = rema + gamma{1,i}*X_real(:,k+1-i);
    end
    X_real(:,k) = diff_X_real - rema;
    %ʵ�ʹ۲�ֵ
    Z_meas(:,k) = 0.1*X_real(1,k) + 0.3*X_real(2,k) + V_noise(1,k);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------����Ӧ�����׿������˲������ܲ���-----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%��ʼֵ���ã���ʼ������Ϊ�㣩

%FCDKF
X_esti_A = zeros(3,N);            %״̬���Ź���ֵ
P_esti_A = cell(1,N);             %����������
X_esti_A(:,1) = [0;0;0];          %״̬���Ƴ�ʼֵ
P_esti_A{1,1} = P_pred_0;         %��ʼ���Ʒ�����

%FEKF
X_esti_B= zeros(3,N);             %״̬���Ź���ֵ
P_esti_B = cell(1,N);             %����������
X_esti_B(:,1) = [0;0;0];          %״̬���Ƴ�ʼֵ
P_esti_B{:,1} = P_pred_0;         %��ʼ���Ʒ�����
%��������δ֪Ϊ��r_test
r_test = 2;
q_test = [0,0,0]';
R_test = 0.3;

r_esti_A = zeros(1,N);            %����������ֵ����ֵ
r_esti_A(1,1) = 2;  
q_esti_A = zeros(3,N);            %ϵͳ������ֵ����ֵ
q_esti_A(:,1) = [1,1,1]';  
R_esti_A = zeros(1,N);            %������������������ֵ                              
R_esti_A(1,1) = 10;    
% Q_esti_A = cell(1,N);           %ϵͳ��������������ֵ                              
% Q_esti_A{1,1} = [0,0,0; 0,0,0; 0,0,0];    

for k=2:1:N
%�������˲�
        %״̬Ԥ��:X_pre
        diff_X_esti = [X_esti_A(2,k-1); ...
                       -0.1*X_esti_A(1,k-1)-X_esti_A(2,k-1)*X_esti_A(3,k-1)+ ...
                       U_input(1,k-1); 0];
            %��������
            rema = [0;0;0];
            if k>L
                for i = 2:1:L+1
                   rema = rema + gamma{1,i}*X_esti_A(:,k+1-i);
                end
            else
                for i = 2:1:k
                    rema = rema + gamma{1,i}*X_esti_A(:,k+1-i);
                end
            end
        X_pre = diff_X_esti - rema + q;     %_esti_A(:,k-1)
        
        %Ԥ�����Э�������:P_pred

            %�����������cholsky�ֽ�
            S_chol = chol(P_esti_A{1,k-1})';

            %��������
            rema_P = [0,0,0;0,0,0;0,0,0];
            if k>L+1
                for i = 2:1:L+2
                    rema_P = rema_P + gamma{1,i}*P_esti_A{1,k+1-i}*gamma{1,i}';
                end
            else
                for i = 2:1:k
                    rema_P = rema_P + gamma{1,i}*P_esti_A{1,k+1-i}*gamma{1,i}';
                end
            end

        temp_f_plus_hs1 = [X_esti_A(2,k-1)+ h_con*S_chol(2,1); ...
                           -0.1*(X_esti_A(1,k-1)+h_con*S_chol(1,1))- ...
                           (X_esti_A(2,k-1)+h_con*S_chol(2,1))*( X_esti_A(3,k-1) ...
                           +h_con*S_chol(3,1))+U_input(1,k-1); 0];
        temp_f_subtraction_hs1 = [X_esti_A(2,k-1)-h_con*S_chol(2,1); ...
                         -0.1*(X_esti_A(1,k-1)-h_con*S_chol(1,1))- ...
                         (X_esti_A(2,k-1)-h_con*S_chol(2,1))* ...
                         (X_esti_A(3,k-1)-h_con*S_chol(3,1))+U_input(1,k-1); 0];   

        temp_f_plus_hs2 = [X_esti_A(2,k-1)+ h_con*S_chol(2,2); ...
                         -0.1*(X_esti_A(1,k-1)+h_con*S_chol(1,2))- ...
                         (X_esti_A(2,k-1)+h_con*S_chol(2,2))* ...
                         (X_esti_A(3,k-1)+h_con*S_chol(3,2))+U_input(1,k-1); 0];
        temp_f_subtraction_hs2 = [X_esti_A(2,k-1)-h_con*S_chol(2,2); ...
                         -0.1*(X_esti_A(1,k-1)-h_con*S_chol(1,2))- ...
                         (X_esti_A(2,k-1)-h_con*S_chol(2,2)) * ...
                         (X_esti_A(3,k-1)-h_con*S_chol(3,2))+U_input(1,k-1); 0]; 

        temp_f_plus_hs3 = [X_esti_A(2,k-1)+h_con*S_chol(2,3); ...
                         -0.1*(X_esti_A(1,k-1)+h_con*S_chol(1,3))- ...
                         (X_esti_A(2,k-1)+h_con*S_chol(2,3))* ...
                         (X_esti_A(3,k-1)+h_con*S_chol(3,3))+U_input(1,k-1); 0];
        temp_f_subtraction_hs3 = [X_esti_A(2,k-1)-h_con*S_chol(2,3); ...
                         -0.1*(X_esti_A(1,k-1)-h_con*S_chol(1,3))- ...
                         (X_esti_A(2,k-1)-h_con*S_chol(2,3)) * ...
                         (X_esti_A(3,k-1)-h_con*S_chol(3,3))+U_input(1,k-1); 0]; 

        temp_initial = ( 0.5*(temp_f_plus_hs1-temp_f_subtraction_hs1)*S_chol(:,1)'+ ...
                         0.5*(temp_f_plus_hs2-temp_f_subtraction_hs2)*S_chol(:,2)'+ ...
                         0.5*(temp_f_plus_hs3-temp_f_subtraction_hs3)*S_chol(:,3)' ...
                         )*gamma{1,2}'; 

        %F = [0,1,0; -0.1,-X_esti(3,k-1),-X_esti(2,k-1); 0,0,0];
        %P_pred = (F-gamma{1,2})*P_esti{1,k-1}*(F-gamma{1,2})'+ Q + rema_P;
        P_pred = 1/(4*h_con^2)*( ...
                  (temp_f_plus_hs1 - temp_f_subtraction_hs1)*(temp_f_plus_hs1 - temp_f_subtraction_hs1)'  ...
                 +(temp_f_plus_hs2 - temp_f_subtraction_hs2)*(temp_f_plus_hs2 - temp_f_subtraction_hs2)'  ...
                 +(temp_f_plus_hs3 - temp_f_subtraction_hs3)*(temp_f_plus_hs3 - temp_f_subtraction_hs3)'  ...
                ) + ...
                - 1/h_con * temp_initial - 1/h_con * temp_initial' + ...
                + Q+ rema_P; %_esti_A{1,k-1} 

        %��״̬Ԥ�����Э����������cholsky�ֽ�
        S_chol_pred = chol(P_pred)';

        %����ֵ����  Z_esti ---- Z_k|k-1
        Z_esti = 0.1*X_pre(1,1) + 0.3*X_pre(2,1) + r_esti_A(1,k-1);%

        %����Ԥ�����Э�������:P_zpred ---- P_z_k|k-1
        P_zpred = 1/(4*h_con^2) * ( ... 
                  ( 0.1* (X_pre(1,1)+h_con*S_chol_pred(1,1)) + 0.3*(X_pre(2,1)+h_con*S_chol_pred(2,1)) ...
                   - ...
                   0.1* (X_pre(1,1)-h_con*S_chol_pred(1,1)) - 0.3*(X_pre(2,1)-h_con*S_chol_pred(2,1)) ...
                   )^2 + ...
                  ...
                  ( 0.1* (X_pre(1,1)+h_con*S_chol_pred(1,2)) + 0.3*(X_pre(2,1)+h_con*S_chol_pred(2,2)) ...
                   - ...
                   0.1* (X_pre(1,1)-h_con*S_chol_pred(1,2)) - 0.3*(X_pre(2,1)-h_con*S_chol_pred(2,2)) ...
                   )^2 + ...
                   ...
                  ( 0.1* (X_pre(1,1)+h_con*S_chol_pred(1,3)) + 0.3*(X_pre(2,1)+h_con*S_chol_pred(2,3)) ...
                   - ...
                   0.1* (X_pre(1,1)-h_con*S_chol_pred(1,3)) - 0.3*(X_pre(2,1)-h_con*S_chol_pred(2,3)) ...
                   )^2 ...
                  ) ...
                  + R;%_esti_A(1,k-1)

        %Ԥ����Э�������:P_xzpred
        Y_pre =  X_pre; % inv(S_chol_pred) * 
        P_xzpred = 1/(2*h_con) * ( ... 
                  ( 0.1* (X_pre(1,1)+h_con*S_chol_pred(1,1)) + 0.3*(X_pre(2,1)+h_con*S_chol_pred(2,1)) ...
                   - ...
                   0.1* (X_pre(1,1)-h_con*S_chol_pred(1,1)) - 0.3*(X_pre(2,1)-h_con*S_chol_pred(2,1)) ...
                   ) * S_chol_pred(:,1) + ...
                  ...
                  ( 0.1* (X_pre(1,1)+h_con*S_chol_pred(1,2)) + 0.3*(X_pre(2,1)+h_con*S_chol_pred(2,2)) ...
                   - ...
                   0.1* (X_pre(1,1)-h_con*S_chol_pred(1,2)) - 0.3*(X_pre(2,1)-h_con*S_chol_pred(2,2)) ...
                   ) * S_chol_pred(:,2) + ...
                   ...
                  ( 0.1* (X_pre(1,1)+h_con*S_chol_pred(1,3)) + 0.3*(X_pre(2,1)+h_con*S_chol_pred(2,3)) ...
                   - ...
                   0.1* (X_pre(1,1)-h_con*S_chol_pred(1,3)) - 0.3*(X_pre(2,1)-h_con*S_chol_pred(2,3)) ...
                   ) * S_chol_pred(:,3)  ...
                   );

        %���㿨��������:Kk(2*1)
        Kk = P_xzpred/P_zpred;

        %״̬����
        X_esti_A(:,k) = X_pre + Kk*( Z_meas(1,k) - Z_esti );

        %������ֵ����
        r_esti_A(k) = ( (k-1)*r_esti_A(k-1) + Z_meas(1,k) - (0.1*X_pre(1,1) + 0.3*X_pre(2,1)) )/k;
        %q_esti_A(:,k) = ( (k-1)*q_esti_A(:,k-1) + X_esti_A(:,k) - (diff_X_esti - rema) )/k;
        
        %���Ʒ���������
        P_esti_A{1,k} = P_pred - Kk*P_zpred*Kk';

        %���������������
         %R_esti_A(1,k) = ( (k-1)*R_esti_A(1,k-1) - P_zpred + R_esti_A(1,k-1) +  ...
         %            ( Z_meas(:,k) - Z_esti )^2 )/k
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------��չ�����׿������˲������ܲ���-----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=2:1:N
  %�������˲�
        %״̬Ԥ��:X_pre
        diff_X_esti = [X_esti_B(2,k-1); -0.1*X_esti_B(1,k-1)- ...
                       X_esti_B(2,k-1)*X_esti_B(3,k-1)+U_input(1,k-1); 0];
            %��������
            rema = [0;0;0];
            if k>L
                for i = 2:1:L+1
                   rema = rema + gamma{1,i}*X_esti_B(:,k+1-i);
                end
            else
                for i = 2:1:k
                    rema = rema + gamma{1,i}*X_esti_B(:,k+1-i);
                end
            end
        X_pre = diff_X_esti - rema + q_test;     %һ��״̬Ԥ��
        %Ԥ�����Э�������:P_pred
            %��������
            rema_P = [0,0,0;0,0,0;0,0,0];
            if k>L+1
                for i = 3:1:L+2
                    rema_P = rema_P + gamma{1,i}*P_esti_B{1,k+1-i}*gamma{1,i}';
                end
            else
                for i = 3:1:k
                    rema_P = rema_P + gamma{1,i}*P_esti_B{1,k+1-i}*gamma{1,i}';
                end
            end
            F = [0,1,0; -0.1,-X_esti_B(3,k-1),-X_esti_B(2,k-1); 0,0,0];
        P_pred = (F-gamma{1,2})*P_esti_B{1,k-1}*(F-gamma{1,2})'+ Q + rema_P;
        %���㿨��������:Kk(2*1)
        H = [0.1,0.3,0];
        Kk = P_pred*H'/(H*P_pred*H' + R_test);
        %״̬����
        X_esti_B(:,k) = X_pre + Kk*( Z_meas(k) - 0.1*X_pre(1,1) - ...
                      0.3*X_pre(2,1) - r_test );
        %���Ʒ���������
        P_esti_B{1,k} = (I-Kk*H)*P_pred;
end

%״̬����ͼ
k = 1:1:N;
figure;
plot(k,X_real(3,:),'b',k,X_esti_A(3,:),'--r',k,X_esti_B(3,:),':g','linewidth',LineWidth);
legend('ʵ��״̬3','FCDKF����״̬3','FEKF����״̬3','Location','best');

figure;
plot(k,X_real(2,:),'b',k,X_esti_A(2,:),'--r',k,X_esti_B(2,:),':g','linewidth',LineWidth);
legend('ʵ��״̬2','FCDKF����״̬2','FEKF����״̬2','Location','best');

figure;
plot(k,X_real(1,:),'b',k,X_esti_A(1,:),'--r',k,X_esti_B(1,:),':g','linewidth',LineWidth);
legend('ʵ��״̬1','FCDKF����״̬1','FEKF����״̬1','Location','best');

%save FEKF_AFCDKF_MIMO_COMPARISON_DATA X_real X_esti_A X_esti_B r_esti_A

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ͼ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %������������ͼ
% k = 1:1:N;
% figure;
% plot(k,U_input,'b',k,Z_meas,':r','linewidth',LineWidth);
% legend('����','�������','Location','best');
% 
% %״̬����ͼ
% figure;
% plot(k,X_real(1,:),'b',k,X_esti_A(1,:),':r','linewidth',LineWidth);
% legend('ʵ��״̬1','����״̬1','Location','best');
% 
% %״̬����ͼ
% figure;
% plot(k,X_real(2,:),'b',k,X_esti_A(2,:),':r','linewidth',LineWidth);
% legend('ʵ��״̬2','����״̬2','Location','best');
% 
% %״̬����ͼ
% figure;
% plot(k,X_real(3,:),'b',k,X_esti_A(3,:),':r','linewidth',LineWidth);
% legend('ʵ��״̬3','����״̬3','Location','best');
% 
% %���ͼ
%������
% err_state = zeros(2,N);
figure;
err_state = X_real - X_esti_A;
%plot(k,q_esti_A(1,k),'b',k,q_esti_A(1,k),'--r',k,q_esti_A(1,k),'.g','linewidth',LineWidth);
plot(k,r_esti_A(1,k),'b','linewidth',LineWidth);
legend('���Ʋ���������ֵ','Location','best');
line([0,N],[r,r]);
%line([0,N],[R,R]);,k,R_esti_A(1,k),'r'






