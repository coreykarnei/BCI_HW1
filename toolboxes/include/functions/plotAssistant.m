p%% This code is very helpful in specifying the exact width of the line and
%% axis. 

% Author : Mohit Kumar GOEL, 7 DEC 2009
% 

% set(findobj(gca,'Type','line'),'LineWidth',2)
% set(get(gca, 'Title'),'fontsize',16,'fontweight','b');
% set(gca, 'LineWidth',2,'FontSize',12,'fontweight','b' )
% set(get(gca, 'XLabel'), 'LineWidth',2,'FontSize',12,'fontweight','b')
% set(get(gca, 'YLabel'), 'LineWidth',2,'FontSize',12,'fontweight','b')
% grid on


%%
ha = get(gca, 'Children');
hf = get(gcf, 'Children');
% Setting the parameters of the axis and the legend size
for i = 1:length(ha)
    set(ha(i), 'LineWidth', 2); % making plot lines thicker
end
for i = 1:length(hf)
    set(hf(i), 'FontSize',12,'fontweight','b') %Setting the size of legend fonts and the axis size
end
% set(findobj(gca,'Type','line'),'LineWidth',2)
set(get(gca, 'Title'),'fontsize',16,'fontweight','b');
set(gca, 'LineWidth',2,'FontSize',12,'fontweight','b' )
% legend('Location','NorthEast')   % this command is not very trustworthy
set(get(gca, 'XLabel'), 'LineWidth',2,'FontSize',12,'fontweight','b')
set(get(gca, 'YLabel'), 'LineWidth',2,'FontSize',12,'fontweight','b')
grid on

