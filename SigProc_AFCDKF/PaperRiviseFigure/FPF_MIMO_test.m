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
I = eye(3,3);                %���ɵ�λ��
%I(3,3) = 0;

%����
q = [0 0 0]';                      %ϵͳ������ֵ
r = 1;                             %����������ֵ
Q = [0.3,0,0; 0,0.3,0; 0,0,0.001]; %ϵͳ�����������
R = 0.3;                           %���������������
W_noise = Q*randn(3,N) + q;        %ϵͳ����
V_noise = R*randn(1,N) + r;        %��������

x = zeros(3,tf); % ϵͳ״̬��ʵֵ ��ʼֵ0
P_pred_0 = [100,0,0; 0,100,0;0,0,100];%��ʼԤ�ⷽ����
x_0  = [0; 0; 0.2];                   %��ʼ״̬
x(:,1) = x_0;                    %��ʵ״̬��ʼֵ
y = zeros(1,tf); % ϵͳ״̬��ʵֵ ��ʼֵ0
y(1,1) = 0.1*x(1,1) + 0.3*x(2,1) +  sqrt(R) * randn + r;

P = cell(1,tf);        % ��������
P{1,1} = P_pred_0;     %��ʼ���Ʒ�����
xhatPart = zeros(3,tf);%״̬����ֵ

%��¼����ʱ���Լ����ķ���
FPF_SISO_ERROR_TIME = zeros(3,50);

xpart = zeros(3,N,tf);
for i = 1 : N
    xpart(:,i,1) = x(:,1) + sqrt(P{1,1}) * randn + q;%��ʼ״̬����x=0��ֵ������Ϊsqrt(P)�ĸ�˹�ֲ�
end
% xArr = [x];
% yArr = [];
% xhatArr = [x];
% PArr = [P];
%xhatPartArr = [xhatPart]; %

%����alpha�״ζ�Ӧ��GL����ϵ�� binomial coefficient
bino_fir = zeros(1,tf);       %΢�ֽ״�Ϊ0.7ʱGL�����µ�ϵ��
alpha = 1;
bino_fir(1,1) = 0.7;
for i = 2:1:tf
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);
end

bino_sec = zeros(1,tf);          %΢�ֽ״�Ϊ1.2ʱGL�����µ�ϵ��
alpha = 1.2;
bino_sec(1,1) = 1;
for i = 2:1:tf
    bino_sec(1,i) = (1-(alpha+1)/(i-1))*bino_sec(1,i-1);  
end

bino_thi = zeros(1,tf);          %΢�ֽ״�Ϊ1ʱGL�����µ�ϵ��
alpha = 1;
bino_thi(1,1) = 1;
for i = 2:1:tf
    bino_thi(1,i) = (1-(alpha+1)/(i-1))*bino_thi(1,i-1);  
end

%����GL΢�ֶ�����ϵ������
gamma = cell(1,tf);
temp_matrx = zeros(3,3);
for i = 1:1:tf 
    temp_matrx(1,1) = bino_fir(1,i);
    temp_matrx(2,2) = bino_sec(1,i);
    temp_matrx(3,3) = bino_thi(1,i);
    gamma{1,i} = temp_matrx;
end

%ϵͳ��������:ǰһ������Ϊ1����һ������Ϊ
U_input = ones(1,tf);            
for i = N/2+1:1:tf
    U_input(:,i) = -2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diff_X_real    ��ʾkʱ��״̬��΢��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff_X_real = [0;0;0];
xpartminus = zeros(3,N);

for k = 2 : tf
    
    diff_X_real = [x(2,k-1); ...
                  -0.1*x(1,k-1)-x(2,k-1)*x(3,k-1)+ ...
                  U_input(1,k-1); 0] + W_noise(:,k-1);
    %3*sin(2*x(1,k-1)) -x(1,k-1) + W_noise(1,k-1);
    rema = [0;0;0];
    for i = 2:1:k
        rema = rema + gamma{1,i}*x(1,k+1-i);
    end
    x(1,k) = diff_X_real - rema;
    %kʱ����ʵֵ
    y(1,k) = 0.1*x(1,k) + 0.3*x(2,k) + V_noise(1,k);  %kʱ�̹۲�ֵ

 %% ����N������
 for i = 1 : N
     %�������N������
     xpartminus(:,i) = [xpart(2,i,k-1); ...
                        -0.1*xpart(1,i,k-1)-xpart(2,i,k-1)*xpart(3,i,k-1)+ ...
                        U_input(1,k-1); 0] + sqrt(Q) * randn;
     %3*sin(2*xpart(i,k-1)) - xpart(i,k-1) + sqrt(Q) * randn;
     temp = 0;
         for j = 2 : 1 : k
            temp = temp + gamma{1,i}*xpart(:,i,k+1-j);
         end
     xpartminus(:,i) = xpartminus(:,i) - temp;
     ypart = xpartminus(:,i);      %ÿ�����Ӷ�Ӧ�۲�ֵ
     vhat = y(1,k) - ypart;        %����ʵ�۲�֮�����Ȼ
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







