
% ����������ѯ��MATLAB����λ��Ƽ�ͷ��Ϊ�˾���һ�������Ƽ�ͷ��ר������ӣ�������һ
% ���ܽᡣ��������Matlab���Դ��Ļ��Ƽ�ͷ�ĺ������о������ܣ��������������Ա������
% ɫɫ�ĺ�������Ͼ�˼��ˡ�
һ����ά��ͷ
    1������annotation�������ƶ�ά��ͷannotation���������ڵ�ǰͼ�δ��ڽ���ע�Ͷ���
    ��annotation���󣩣����ĵ��ø�ʽ���£�
(1) annotation(annotation_type)  % ��ָ���Ķ������ͣ�ʹ��Ĭ������ֵ����ע�Ͷ���
(2) annotation('line',x,y)       % ������(x(1), y(1))��(x(2), y(2))����ע�Ͷ���
(3) annotation('arrow',x,y)      % ������(x(1), y(1))��(x(2), y(2))�ļ�ͷע�Ͷ���
(4) annotation('doublearrow',x,y)% ������(x(1), y(1))��(x(2), y(2))��˫��ͷע�Ͷ���
(5) annotation('textarrow',x,y)  % ������(x(1),y(1))��(x(2),y(2))�Ĵ��ı���ļ�ͷע�Ͷ���
(6) annotation('textbox',[x y w h])  % �����ı���ע�Ͷ������½�����(x,y)����w����h.
(7) annotation('ellipse',[x y w h])  % ������Բ��ע�Ͷ���
(8) annotation('rectangle',[x y w h])% ��������ע�Ͷ���
(9) annotation(figure_handle,��)     % �ھ��ֵΪfigure_handle��ͼ�δ��ڽ���ע�Ͷ���
(10)annotation(��,'PropertyName',PropertyValue,��)  % ����������ע�Ͷ�������ԡ�
(11)anno_obj_handle = annotation(��)  % ����ע�Ͷ���ľ��ֵ��
% ע�⣺annotation����ĸ�������figure���������ᵽ������x��y�Ǳ�׼�������꣬��
% ����ͼ�δ��ڣ�figure�������½�Ϊ(0,  0)�����Ͻ�Ϊ(1,  1)�����w�͸߶�hҲ����
% ��׼���ģ���ȡֵ��[0,  1]֮�䡣

 P = [3 1; 1 4]; 
 r = 5;
 [V, D] = eig(P);     % ������ֵ������Բ��Ϊ��׼����
 a = sqrt(r/D(1));     % ��Բ������
 b = sqrt(r/D(4));    % ��Բ�̰���
 t = linspace(0, 2*pi, 60);    % �ȼ������һ����0��2pi�İ���60��Ԫ�ص�����
 xy = V*[a*cos(t); b*sin(t)];    % ������Բ�ļ����귽�̼�����Բ�ϵ������
 plot(xy(1,:),xy(2,:), 'k', 'linewidth', 3);    % ������Բ���ߣ��߿�Ϊ3����ɫΪ��ɫ
% �ڵ�ǰͼ�δ��ڼ������ͷ���ı���ע��
 h = annotation('textarrow',[0.606 0.65],[0.55 0.65]);
% �����ı���ע������ʾ���ַ����������ֺ�Ϊ15
 set(h, 'string','3x^2+2xy+4y^2 = 5', 'fontsize', 15);
 annotation('doublearrow',[0.2 0.8],[0.85 0.85],...
        'LineStyle','-','color',[1 0 0],'HeadStyle','cback3');
    
    
% ���Ƶ����ǣ���������ǵ�λ��
>> cla reset;
>> load topo;
>> [x y z] = sphere(45);
>> s = surface(x,y,z,'FaceColor','texturemap','CData',topo);
>> colormap(topomap1);
% Brighten the colormap for better annotation visibility:
>> brighten(.6)
% Create and arrange the camera and lighting for better visibility:
>> campos([1.3239  -14.4250  9.4954]);
>> camlight;
>> lighting gouraud;
>> axis off vis3d;
% Set the x- and y-coordinates of the textarrow object:
>> x = [0.7698 0.5851];
>> y = [0.3593 0.5492];
% Create the textarrow object: 
>> txtar =  annotation('textarrow',x,y,'String','We are here.','FontSize',14);
  

2������quiver�������Ƽ�ͷ
quiver�����ĵ��ø�ʽ���£�
quiver(x,y,u,v)
quiver(u,v)
quiver(...,scale)
quiver(...,LineSpec)
quiver(...,LineSpec,'filled')
quiver(axes_handle,...)
h = quiver(...)

����3�������������ߣ�������ͼ�Ρ� 
    


