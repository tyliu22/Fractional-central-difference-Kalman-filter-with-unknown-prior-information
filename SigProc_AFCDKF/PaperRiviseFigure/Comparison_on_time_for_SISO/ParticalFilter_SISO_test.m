%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   �����������˲����渴��
%   ���ģ�     fractional order PF
%   Ŀ�ģ������������˲��㷨����
%         ��ϵͳ������ֵ���й���
%         ����ʵ��:    D^{0.7} x_k = 3*sin(2*x_{k-1}) -x_{k-1} + w_k
%                              y_k = x_k + v_k
%   ������ϺõĶ�״̬���й��ƣ���ֵϵͳ������ֵ����
%
%   ��ע�������������˲����㷨����
%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% clc;
% clear all;

clc;
clear;
tempp = 1;

tic

h_con = sqrt(3);
tf = 100; % ����ʱ��
N = 50;  % ���Ӹ���

%ϵͳ��������
% A = [0,1; -0.1,-0.2];      %ϵͳ����
% B = [0; 1];                %
% C = [0.1,0.3];             %
I = eye(1,1);                %���ɵ�λ��
%I(3,3) = 0;

%����
q = 0;                      %ϵͳ������ֵ
r = 0;                      %����������ֵ
Q = 0.81;                   %ϵͳ�����������
R = 0.25;                   %���������������
W_noise = sqrt(Q)*randn(1,tf) + q;  %ϵͳ����
V_noise = sqrt(R)*randn(1,tf) + r;  %��������
x = zeros(1,tf); % ϵͳ״̬��ʵֵ ��ʼֵ0
y = zeros(1,tf); % ϵͳ״̬��ʵֵ ��ʼֵ0
y(1,1) = x(1,1) + sqrt(R) * randn + r;

P = zeros(1,tf); % ��������
P(1,1) = 2;      % ��ʼ�����ֲ��ķ���
xhatPart = zeros(1,tf);%״̬����ֵ

%��¼����ʱ���Լ����ķ���
FPF_SISO_ERROR_TIME = zeros(3,50);

xpart = zeros(N,tf);
for i = 1 : N
    xpart(i,1) = x(1,1) + sqrt(P(1,1)) * randn + q;%��ʼ״̬����x=0��ֵ������Ϊsqrt(P)�ĸ�˹�ֲ�
end
% xArr = [x];
% yArr = [];
% xhatArr = [x];
% PArr = [P];
%xhatPartArr = [xhatPart]; %

%����alpha�״ζ�Ӧ��GL����ϵ�� binomial coefficient
bino_fir = zeros(1,tf);       %΢�ֽ״�Ϊ0.7ʱGL�����µ�ϵ��
alpha = 0.7;
bino_fir(1,1) = 1;
for i = 2:1:tf
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diff_X_real    ��ʾkʱ��״̬��΢��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff_X_real = 0;

for k = 2 : tf

    diff_X_real = 3*sin(2*x(1,k-1)) -x(1,k-1) + W_noise(1,k-1);
    rema = 0;
    for i = 2:1:k
        rema = rema + bino_fir(1,i)*x(1,k+1-i);
    end
    x(1,k) = diff_X_real - rema;
    %kʱ����ʵֵ
    y(1,k) = x(1,k) + V_noise(1,k);  %kʱ�̹۲�ֵ

 %% ����N������
 for i = 1 : N
     %�������N������
     xpartminus(i) = 3*sin(2*xpart(i,k-1)) - xpart(i,k-1) + sqrt(Q) * randn;
     temp = 0;
         for j = 2 : 1 : k
            temp = temp + bino_fir(1,j)*xpart(i,k+1-j);
         end
     xpartminus(i) = xpartminus(i) - temp;
     ypart = xpartminus(i);      %ÿ�����Ӷ�Ӧ�۲�ֵ
     vhat = y(1,k) - ypart;             %����ʵ�۲�֮�����Ȼ
     q(i) = (1 / sqrt(R) / sqrt(2*pi)) * exp(-vhat^2 / 2 / R);
     %ÿ�����ӵ���Ȼ�����ƶ�
 end

 %%
%Ȩֵ��һ��
qsum = sum(q);
for i = 1 : N
    q(i) = q(i) / qsum; %��һ�����Ȩֵ q
end

%%
 %����Ȩֵ���²���
  for i = 1 : N 
      u = rand;
      qtempsum = 0; 
      for j = 1 : N
          qtempsum = qtempsum + q(j); 
          if qtempsum >= u 
              xpart(i,k) = xpartminus(j);
              break;
          else
              xpart(i,k) = xpart(i,k-1);
          end 
      end
  end
xhatPart(1,k) = mean(xpart(:,k));

%%
%����״̬����ֵ��ΪN�����ӵ�ƽ��ֵ�����ﾭ�����²�����������ӵ�Ȩֵ��ͬ
% xArr = [xArr x];   
% yArr = [yArr y];  
% % xhatArr = [xhatArr xhat]; 
% PArr = [PArr P]; 
% xhatPartArr = [xhatPartArr xhatPart];

end

 FPF_SISO_TIME(1,tempp)    = toc
 FPF_ERROR_norm1(1,tempp)  = norm((x - xhatPart),1);
 FPF_ERROR_norm2(1,tempp)  = norm((x - xhatPart),2);
 tempp = tempp + 1;
 
%  FPF_SISO_ERROR_TIME(1,:)  = FPF_ERROR_norm1(1,:);
%  FPF_SISO_ERROR_TIME(2,:)  = FPF_ERROR_norm2(1,:);
%  FPF_SISO_ERROR_TIME(3,:)  = FPF_SISO_TIME(1,:);
%  FPF_ERROR_norm1_average  = sum(FPF_SISO_ERROR_TIME(1,:))/50
%  FPF_ERROR_norm2_average  = sum(FPF_SISO_ERROR_TIME(2,:))/50
%  FPF_SISO_TIME_average    = sum(FPF_SISO_ERROR_TIME(3,:))/50
%  save FPF_SISO_ERROE_TIME FPF_SISO_ERROR_TIME FPF_ERROR_norm1_average ...
%       FPF_ERROR_norm2_average FPF_SISO_TIME_average







