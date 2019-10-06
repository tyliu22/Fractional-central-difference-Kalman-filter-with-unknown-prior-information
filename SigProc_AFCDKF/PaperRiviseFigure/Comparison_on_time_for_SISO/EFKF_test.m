%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   �����׿������˲������渴��
%   ���ģ�     fractional order CDKF
%   Ŀ�ģ�һ�׷��������Ĳ�ֿ������˲����㷨����
%         ��ϵͳ������ֵ���й���
%         ����ʵ��:    D^{0.7} x_k = 3*sin(2*x_{k-1}) -x_{k-1} + w_k
%                              y_k = x_k + v_k
%   ������ϺõĶ�״̬���й��ƣ���ֵϵͳ������ֵ����
%
%   ��ע����������չ�������˲������㷨����
%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  clc
%  clear
%  tempp = 1;

 tic

%���沽��
N = 50;

%״̬������ʼ��
X_real = zeros(1,N);         %��ʵ״̬
X_esti = zeros(1,N);         %״̬���Ź���ֵ
P_esti = zeros(1,N);         %����������
Z_meas = zeros(1,N);         %ʵ�ʹ۲�ֵ

%ϵͳ��������
% A = [0,1; -0.1,-0.2];      %ϵͳ����
% B = [0; 1];                %
% C = [0.1,0.3];             %
I = eye(1,1);                %���ɵ�λ��
%I(3,3) = 0;

%����
q = 0;                     %ϵͳ������ֵ
r = 0;                     %����������ֵ
Q = 0.81;                    %ϵͳ�����������
R = 0.25;                    %���������������
W_noise = sqrt(Q)*randn(1,N) + q;  %ϵͳ����
V_noise = sqrt(R)*randn(1,N) + r;  %��������

%��ʼֵ����
P_pred_0 = 100;              %��ʼԤ�ⷽ����
P_esti(1,1) = P_pred_0;      %��ʼ���Ʒ�����
x_0  = 0;                    %��ʼ״̬
X_real(1,1) = x_0;           %��ʵ״̬��ʼֵ
X_esti(1,1) = 0;             %״̬���Ƴ�ʼֵ
Z_meas(1,1) = V_noise(1,1);  %�������ݳ�ʼֵ

FEKF_SISO_ERROR_TIME = zeros(3,50);%һ����50�η�������

%����alpha�״ζ�Ӧ��GL����ϵ�� binomial coefficient 
bino_fir = zeros(1,N);       %΢�ֽ״�Ϊ0.7ʱGL�����µ�ϵ��
alpha = 0.7;
bino_fir(1,1) = 1;
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
        diff_X_esti = 3*sin(2*X_esti(1,k-1)) - X_esti(1,k-1) + q;
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
        X_pre = diff_X_esti - rema;     %һ��״̬Ԥ��
        %Ԥ�����Э�������:P_pred
            %��������
            rema_P = 0;
            if k>L+1
                for i = 3:1:L+2
                    rema_P = rema_P + bino_fir(1,i)*P_esti(1,k+1-i)*bino_fir(1,i)';
                end
            else
                for i = 3:1:k
                    rema_P = rema_P + bino_fir(1,i)*P_esti(1,k+1-i)*bino_fir(1,i)';
                end
            end
            F = 6*cos(2*X_esti(1,k-1)) - 1;
        P_pred = (F-bino_fir(1,2))*P_esti(1,k-1)*(F-bino_fir(1,2))'+ Q + rema_P;
        
        %���㿨��������:Kk(2*1)
        H = 1;
        Kk = P_pred*H'/(H*P_pred*H' + R);
        
        %״̬����
        X_esti(1,k) = X_pre + Kk*( Z_meas(1,k) - X_pre -r);
        
        %���Ʒ���������
        P_esti(1,k) = (I-Kk*H)*P_pred;

end

 FEKF_SISO_TIME(1,tempp)    = toc
 FEKF_ERROR_norm1(1,tempp)  = norm((X_real - X_esti),1)
 FEKF_ERROR_norm2(1,tempp)  = norm((X_real - X_esti),2)
 tempp = tempp + 1;
 
%  FEKF_SISO_ERROR_TIME(1,:) = FEKF_ERROR_norm1(1,:);
%  FEKF_SISO_ERROR_TIME(2,:) = FEKF_ERROR_norm2(1,:);
%  FEKF_SISO_ERROR_TIME(3,:) = FEKF_SISO_TIME(1,:);
%  FEKF_ERROR_norm1_average  = sum(FEKF_SISO_ERROR_TIME(1,:))/50
%  FEKF_ERROR_norm2_average  = sum(FEKF_SISO_ERROR_TIME(2,:))/50
%  FEKF_SISO_TIME_average    = sum(FEKF_SISO_ERROR_TIME(3,:))/50
%  save FEKF_SISO_ERROE_TIME1 FEKF_SISO_ERROR_TIME FEKF_ERROR_norm1_average ...
%       FEKF_ERROR_norm2_average FEKF_SISO_TIME_average










