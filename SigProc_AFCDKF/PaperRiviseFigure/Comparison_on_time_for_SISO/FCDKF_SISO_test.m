%
%   �����׿������˲������渴��
%   ���ģ�     fractional order CDKF
%   Ŀ�ģ�һ�׷��������Ĳ�ֿ������˲����㷨����
%         ��ϵͳ������ֵ���й���
%         ����ʵ��:    D^{0.7} x_k = 3*sin(2*x_{k-1}) -x_{k-1} + w_k
%                              y_k = x_k + v_k
%   �����
%
%   ��ע�����������Ĳ�ֿ������˲����㷨
%         *********����ʱ�����*********
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%  clc
%  clear
%  tempp = 1;

tic
 
%���沽��
N = 100;
h_con = sqrt(2.0);

%״̬������ʼ��
X_real = zeros(1,N);         %��ʵ״̬
X_esti = zeros(1,N);         %״̬���Ź���ֵ
P_esti = zeros(1,N);         %����������
Z_meas = zeros(1,N);         %ʵ�ʹ۲�ֵ

%ϵͳ��������
I = eye(1,1);                %���ɵ�λ��
 
%����
q = 0;                      %ϵͳ������ֵ
r = 0;                      %����������ֵ
Q = 0.81;                   %ϵͳ�����������
R = 0.25;                   %���������������
W_noise = sqrt(Q)*randn(1,N) + q;  %ϵͳ����
V_noise = sqrt(R)*randn(1,N) + r;  %��������

%��ʼֵ���ã���ʼ������Ϊ�㣩
P_pred_0 = 100;              %��ʼԤ�ⷽ����
P_esti(1,1) = P_pred_0;      %��ʼ���Ʒ�����
x_0  = 0;                    %��ʼ״̬     
X_real(1,1) = x_0;           %��ʵ״̬��ʼֵ
X_esti(1,1) = 0;             %״̬���Ƴ�ʼֵ
Z_meas(1,1) = V_noise(1,1);  %�������ݳ�ʼֵ

%FCDKF_SE = zeros(1,N);
FCDKF_SISO_ERROR_TIME = zeros(3,50);%һ����50�η�������

%����alpha�״ζ�Ӧ��GL����ϵ�� binomial coefficient 
bino_fir = zeros(1,N);       %΢�ֽ״�Ϊ0.7ʱGL�����µ�ϵ��
alpha = 1;
bino_fir(1,1) = 0.7;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);  
end

%GL�����¶̼���ԭ��ĳ���
L = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diff_X_real    ��ʾkʱ��״̬��΢��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff_X_real = 0;

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
  %�������˲�
        %״̬Ԥ��:X_pre
        diff_X_esti = 3*sin(2*X_esti(1,k-1)) - X_esti(1,k-1);
            %��������
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
        X_pre = diff_X_esti - rema + q;         %һ��״̬Ԥ��
        %Ԥ�����Э�������:P_pred

            %�����������cholsky�ֽ�
            S_chol = chol(P_esti(1,k-1))';

            %��������
            rema_P = 0;
            if k>L+1
                for i = 2:1:L+2
                    rema_P = rema_P + bino_fir(1,i)*P_esti(1,k+1-i)*bino_fir(1,i)';
                end
            else
                for i = 2:1:k
                    rema_P = rema_P + bino_fir(1,i)*P_esti(1,k+1-i)*bino_fir(1,i)';
                end
            end

        %��ʱ���� temp_fun : ������ֵ,����Ϊ������
        temp_fun = 3*sin( 2*(X_esti(1,k-1)+h_con*S_chol) ) - (X_esti(1,k-1)+h_con*S_chol) - ...
                   3*sin( 2*(X_esti(1,k-1)-h_con*S_chol) ) + (X_esti(1,k-1)-h_con*S_chol);
        P_pred = 1/(4*h_con^2) * temp_fun^2 + rema_P + ...
                 1/h_con*0.5*temp_fun*S_chol'*(-bino_fir(1,2))' + ...
                 1/h_con*(-bino_fir(1,2))*S_chol*0.5*temp_fun' + Q;

        %����ֵ����  Z_esti ---- Z_k|k-1
        Z_esti = X_pre + r;
        
        %����Ԥ�����Э�������:P_zpred ---- P_z_k|k-1
        P_zpred = P_pred + R;

        %���㿨��������:Kk(2*1)
        Kk = P_pred/P_zpred;
        
        %״̬����
        X_esti(1,k) = X_pre + Kk*( Z_meas(1,k) - Z_esti );

        %���Ʒ���������
        P_esti(1,k) = P_pred - Kk*P_zpred*Kk';  
end

 FCDKF_SISO_TIME(1,tempp)    = toc
 FCDKF_ERROR_norm1(1,tempp)  = norm((X_real - X_esti),1)  
 FCDKF_ERROR_norm2(1,tempp)  = norm((X_real - X_esti),2)
 tempp = tempp + 1;

%  FCDKF_SISO_ERROR_TIME(1,:) = FCDKF_ERROR_norm1(1,:);
%  FCDKF_SISO_ERROR_TIME(2,:) = FCDKF_ERROR_norm2(1,:);
%  FCDKF_SISO_ERROR_TIME(3,:) = FCDKF_SISO_TIME(1,:);
%  FCDKF_ERROR_norm1_average  = sum(FCDKF_SISO_ERROR_TIME(1,:))/50
%  FCDKF_ERROR_norm2_average  = sum(FCDKF_SISO_ERROR_TIME(2,:))/50
%  FCDKF_SISO_TIME_average    = sum(FCDKF_SISO_ERROR_TIME(3,:))/50
%  save FCDKF_SISO_ERROE_TIME1 FCDKF_SISO_ERROR_TIME FCDKF_ERROR_norm1_average ...
%       FCDKF_ERROR_norm2_average FCDKF_SISO_TIME_average






