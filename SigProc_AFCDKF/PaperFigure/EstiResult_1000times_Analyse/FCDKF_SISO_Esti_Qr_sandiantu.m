%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ����Ӧ���������Ĳ�ֿ������˲����������ܷ���
%   ���ģ�     fractional order CDKF
%   Ŀ�ģ�ÿ������100��Q��r�Ĺ������ݣ��������� estiqr.mat
%         ����ʵ��:    D^{0.7} x_k = 3*sin(2*x_{k-1}) -x_{k-1} + w_k
%                              y_k = x_k + v_k
%
%   ��ע����Ҫ���SaveData.mһ��ʹ�ã�һ����Ҫ����1000�ι�������    
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

%���沽��
N = 1000;
LineWidth = 2;
h_con = sqrt(1.3);

q_con = 6;                   %ϵͳ������ֵ
r_con = 5;                   %����������ֵ
Q_con = 5;                  %ϵͳ�����������
R_con = 0.25;                %���������������

q = q_con;                   %ϵͳ������ֵ
r = r_con;                   %����������ֵ
Q = Q_con;                   %ϵͳ�����������
R = R_con;                   %���������������

%����alpha�״ζ�Ӧ��GL����ϵ�� binomial coefficient 
bino_fir = zeros(1,N);       %΢�ֽ״�Ϊ0.7ʱGL�����µ�ϵ��
alpha = 0.7;
bino_fir(1,1) = 1;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);  
end

times = 100;
NumOrder = 1;
Q_point = zeros(1,times);
q_point = zeros(1,times);

for kk = 1:1:times
    
    NumOrder
    
    %״̬������ʼ��
    X_real = zeros(1,N);         %��ʵ״̬
    X_esti = zeros(1,N);         %״̬���Ź���ֵ
    P_xesti = zeros(1,N);        %����������
    Z_meas = zeros(1,N);         %ʵ�ʹ۲�ֵ

    %ϵͳ��������
    % A = [0,1; -0.1,-0.2];      %ϵͳ����
    % B = [0; 1];                %
    % C = [0.1,0.3];             %
    I = eye(1,1);                %���ɵ�λ��
    %I(3,3) = 0;

    %����
    W_noise = sqrt(Q)*randn(1,N) + q;  %ϵͳ����
    V_noise = sqrt(R)*randn(1,N) + r;  %��������

    %��ʼֵ���ã���ʼ������Ϊ�㣩
    q_esti = zeros(1,N);         %ϵͳ������ֵ����ֵ
    r_esti = zeros(1,N);         %����������ֵ����ֵ
    Q_esti = zeros(1,N);         %ϵͳ��������������ֵ
    R_esti = zeros(1,N);         %������������������ֵ

    P_pred_0 = 100;              %��ʼԤ�ⷽ����
    P_xesti(1,1) = P_pred_0;     %��ʼ���Ʒ�����
    x_0  = 0;                    %��ʼ״̬     
    X_real(1,1) = x_0;           %��ʵ״̬��ʼֵ
    X_esti(1,1) = 0;             %״̬���Ƴ�ʼֵ
    Z_meas(1,1) = V_noise(1,1);  %�������ݳ�ʼֵ

    %GL�����¶̼���ԭ��ĳ���
    L = N+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % diff_X_real    ��ʾkʱ��״̬��΢��
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    diff_X_real = 0;

    %*************************************************************************%
    %     q_esti(k) = ( (k-1)*q_esti(k-1) + X_esti(1,k) - X_pre(k) + q_esti(k-1) )/k;
    %     Q_esti(k) = ( (k-1)*Q_esti(k-1) + (X_esti(1,k)-X_pre(k))^2 )/k;
    %     r_esti(k) = ( (k-1)*r_esti(k-1) + Z_meas(1,k) - Z_esti )/k;
    %     R_esti(k) = ( (k-1)*R_esti(k-1) + (Z_meas(1,k) - Z_esti)^2 )/k;
    %*************************************************************************%

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
            temp_fun = 3*sin( 2*(X_esti(1,k-1)+h_con*S_chol) ) - (X_esti(1,k-1)+h_con*S_chol) - ...
                       3*sin( 2*(X_esti(1,k-1)-h_con*S_chol) ) + (X_esti(1,k-1)-h_con*S_chol);
            temp(k) = 1/(4*h_con^2) * temp_fun^2 + rema_P + ...
                      1/h_con*0.5*temp_fun*S_chol'*(-bino_fir(1,2))' + ...
                      1/h_con*(-bino_fir(1,2))*S_chol*0.5*temp_fun';
            P_xpred = temp(k) + Q_esti(k-1);
            %����ֵ����  Z_esti ---- Z_k|k-1
            Z_esti = X_pre + r_esti(k-1);

            %����Ԥ�����Э�������:P_zpred ---- P_z_k|k-1
            P_zpred = P_xpred + R;

            %���㿨��������:Kk(2*1)
            Kk = P_xpred/P_zpred;

            %r_esti(k) = ( (k-1)*r_esti(k-1) + Z_meas(1,k) - X_pre )/k;
            r_esti(k) = ( (k-1)*r_esti(k-1) + Z_meas(1,k) - X_pre )/k;

            %״̬����
            X_esti(1,k) = X_pre + Kk*( Z_meas(1,k) - Z_esti );

            %���Ʒ���������
            P_xesti(1,k) = P_zpred - Kk*P_zpred*Kk';

            Q_esti(k) = ( (k-1)*Q_esti(k-1) + ( X_esti(1,k) - X_pre )^2 - temp(k) + P_xesti(1,k) )/k;
            %q_esti(k) = ( (k-1)*q_esti(k-1) +   X_esti(1,k) - diff_X_esti + rema )/k;

    end

    Q_point(1,NumOrder) = Q_esti(k);
    r_point(1,NumOrder) = r_esti(k);

    NumOrder = NumOrder+1;

end

save estiqr Q_point r_point

% ǰһ���������д˶Σ�֮��ֱ������SaveData.m
% Q_data = Q_point;
% r_data = r_point;
% save EstiData Q_data r_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ͼ��� ��ֵ�������ɢ��ͼ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure;
% scatter(r_point(1,:),Q_point(1,:),'b','linewidth',LineWidth);
% legend('��ֵ�������ɢ��ͼ','Location','best');
% hold on
% plot(r_con,Q_con,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],...
%     'MarkerSize',5, 'Marker','o','LineWidth',1, 'Color',[0 0 0]);
% hold off
% 
% figure;
% scatterhist(r_point,Q_point)
% hold on
% plot(r_con,Q_con,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],...
%     'MarkerSize',5, 'Marker','o','LineWidth',1, 'Color',[0 0 0]);
% hold off
% 
% figure;
% scatterhist(r_point,Q_point,'Kernel','on');
% 
