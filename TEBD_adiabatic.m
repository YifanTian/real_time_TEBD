


data1 = load('Data_01_E-8.txt');
data2 = load('Data_015_E-8.txt');
data3 = load('Data_02_E-8.txt');

% plot(data1(:,1),data1(:,2),'-','LineWidth',1.8)
% xlabel('real time step')
% %legend('<Sz1>')
% hold on
% plot(data1(:,1),data1(:,3),'-','LineWidth',1.8)
% xlabel('real time step')
% legend('<Sz1>','<Sz2>')
% hold on
% grid on
plot(data1(:,1),data1(:,4),'-','LineWidth',1.8)
hold on
plot(data2(:,1),data2(:,4),'-','LineWidth',1.8)
hold on
plot(data3(:,1),data3(:,4),'-','LineWidth',1.8)
%xlabel('real time step')
%legend('SvN')
%plot(data1(:,1),data1(:,5),'-','LineWidth',1.8)
xlabel('real time step')
legend('0.01','0.15','0.2')
title('SvN of various tau')
grid on

plot(data1(:,1),data1(:,6),'-','LineWidth',1.8)
hold on
plot(data2(:,1),data2(:,6),'-','LineWidth',1.8)
hold on
plot(data3(:,1),data3(:,6),'-','LineWidth',1.8)
%xlabel('real time step')
%legend('SvN')
%plot(data1(:,1),data1(:,5),'-','LineWidth',1.8)
xlabel('real time step')
legend('0.01','0.15','0.2')
title('SvN of various tau')
grid on


plot(data1(:,1),data1(:,6),'-','LineWidth',1.8)
title('energy bond SS tau = 0.15')
xlabel('real time step')
grid on

fig = figure('Color','w');
%plotyy
%[ax, s1h1 s1h3] = plotyy(data1(:,1),data1(:,2:3),data1(:,1),data1(:,4),'plot');
[ax, s1h1 s1h3] = plotyy(data1(:,1),data1(:,4),data1(:,1),data1(:,5),'plot');
%sig1 color
sig1col = [0 150 150]/255;
%sig1log color
sig1logcol = [210 30 50]/255;
%style the plot
set(s1h1,'Color',sig1col,'LineWidth',3);
set(s1h3,'Color',sig1logcol,'LineWidth',3);
set(ax(1),'YColor',sig1col);
set(ax(2),'YColor',sig1logcol);

%x-axis and y-axis labels
xlabel('$x$','Interpreter','latex');
set(get(ax(1),'Ylabel'),'String','$\mbox{SvN}$','Interpreter','latex')
set(get(ax(2),'Ylabel'),'String','$\mbox{m}$','Interpreter','latex')
%add legend
%leg = legend('$\mbox{$Sz$}$ ','$\mbox{$Sz$}$ ','$\mbox{$SvN$}$ ', 'Location', 'NorthWest');
leg = legend('$\mbox{$SvN$}$ ','$\mbox{$m$}$ ', 'Location', 'NorthWest');

set(leg,'Interpreter','latex');

grid on
