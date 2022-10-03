
function [w,sphere] = eegica(dataica, protocol);
%% runica

% trialICA = 10;
% dataica = data(:,1:trialLeng(trialICA,4));
% dataica = completedata(:,1:round(0.1*size(completedata,2)));

[w,sphere] = runica(dataica,'pca',16);
icacomp = w * data;

figure(2);
for i = 1:16
    subplot(8,2,i)
    plot(icacomp(i,:));
end

if protocol ~= 'MI'
    return;
end

figure;
load MIelectrode;
colorjet=jet;
% minw=min(min(w));
% maxw=max(max(w));
for i = 1:16
    subplot(4,4,i);
    elecMap = zeros(5,4);
    minw=min(w(i,:));
    maxw=max(w(i,:));
    for j = 1:16
        colorindex=round(63*(w(i,j)-minw)/(maxw-minw))+1;
        plot(electrode.position(j,1),electrode.position(j,2),'o',...
            'MarkerEdgeColor','k','MarkerSize',20,...
            'MarkerFaceColor',colorjet(colorindex,:));
        hold on;
        text(electrode.position(j,1) - 0.04 - 0.08*(electrode.position(j,1)>=10), electrode.position(j,2),num2str(j)); hold on;
    end
    title(['Component ' num2str(i)],'fontsize',15);
    grid on;
end