
% 看到有网友询问MATLAB中如何绘制箭头，为此就做一个“绘制箭头”专题的帖子，算是做一
% 个总结。本帖仅对Matlab中自带的绘制箭头的函数进行举例汇总，对于其他网友自编的形形
% 色色的函数就无暇顾及了。
一、二维箭头
    1．调用annotation函数绘制二维箭头annotation函数用来在当前图形窗口建立注释对象
    （annotation对象），它的调用格式如下：
(1) annotation(annotation_type)  % 以指定的对象类型，使用默认属性值建立注释对象。
(2) annotation('line',x,y)       % 建立从(x(1), y(1))到(x(2), y(2))的线注释对象。
(3) annotation('arrow',x,y)      % 建立从(x(1), y(1))到(x(2), y(2))的箭头注释对象。
(4) annotation('doublearrow',x,y)% 建立从(x(1), y(1))到(x(2), y(2))的双箭头注释对象。
(5) annotation('textarrow',x,y)  % 建立从(x(1),y(1))到(x(2),y(2))的带文本框的箭头注释对象
(6) annotation('textbox',[x y w h])  % 建立文本框注释对象，左下角坐标(x,y)，宽w，高h.
(7) annotation('ellipse',[x y w h])  % 建立椭圆形注释对象。
(8) annotation('rectangle',[x y w h])% 建立矩形注释对象。
(9) annotation(figure_handle,…)     % 在句柄值为figure_handle的图形窗口建立注释对象。
(10)annotation(…,'PropertyName',PropertyValue,…)  % 建立并设置注释对象的属性。
(11)anno_obj_handle = annotation(…)  % 返回注释对象的句柄值。
% 注意：annotation对象的父对象是figure对象，上面提到的坐标x，y是标准化的坐标，即
% 整个图形窗口（figure对象）左下角为(0,  0)，右上角为(1,  1)。宽度w和高度h也都是
% 标准化的，其取值在[0,  1]之间。

 P = [3 1; 1 4]; 
 r = 5;
 [V, D] = eig(P);     % 求特征值，将椭圆化为标准方程
 a = sqrt(r/D(1));     % 椭圆长半轴
 b = sqrt(r/D(4));    % 椭圆短半轴
 t = linspace(0, 2*pi, 60);    % 等间隔产生一个从0到2pi的包含60个元素的向量
 xy = V*[a*cos(t); b*sin(t)];    % 根据椭圆的极坐标方程计算椭圆上点的坐标
 plot(xy(1,:),xy(2,:), 'k', 'linewidth', 3);    % 绘制椭圆曲线，线宽为3，颜色为黑色
% 在当前图形窗口加入带箭头的文本标注框
 h = annotation('textarrow',[0.606 0.65],[0.55 0.65]);
% 设置文本标注框中显示的字符串，并设字号为15
 set(h, 'string','3x^2+2xy+4y^2 = 5', 'fontsize', 15);
 annotation('doublearrow',[0.2 0.8],[0.85 0.85],...
        'LineStyle','-','color',[1 0 0],'HeadStyle','cback3');
    
    
% 绘制地球仪，并标出我们的位置
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
  

2．调用quiver函数绘制箭头
quiver函数的调用格式如下：
quiver(x,y,u,v)
quiver(u,v)
quiver(...,scale)
quiver(...,LineSpec)
quiver(...,LineSpec,'filled')
quiver(axes_handle,...)
h = quiver(...)

【例3】绘制正弦曲线，并修饰图形。 
    


