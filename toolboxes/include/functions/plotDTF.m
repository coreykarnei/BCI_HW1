function plotDTF(DTF,fre,dispmodes);

h = gcf;
load options;
load model;

options.dispmodes = dispmodes;

model.electrodes.labels = DTF.labels;
model.electrodes.locations = DTF.locations;
frequency = DTF.frequency;
%%

% get current frequency and matrix
currentf = frequency(1)+fre-1;

currentdtfmatrix = DTF.matrix;%(:,:,fre);
options.srate = 512;  % ??

% enable_specific(true);			% ??
% enable_band(false);			% ??

% update options
options.dtf = 1;

valmin = min(min(currentdtfmatrix));
valmax = max(max(currentdtfmatrix));

% valmin = 0.001;
% valmax = 0.5;

if isnan(valmin)
    valmin = 0.001;
end
if isnan(valmax)
    valmax = 0.5;
end

options.minmax = [valmin, valmax];

options.displimits = [0 valmax];

options.dispmodes = 'out2in';

%% display DTF Graph
hold on;

if options.isskin % draw the skin model
    patch('SpecularStrength',0.2,'DiffuseStrength',0.8,'AmbientStrength',0.8,...
        'FaceLighting','phong',...
        'Vertices',model.skin.italyskin.Vertices,...
        'LineStyle','none',...
        'Faces',model.skin.italyskin.Faces,...
        'FaceColor','interp',...
        'EdgeColor','none',...
        'FaceVertexCData',model.skin.italyskin.FaceVertexCData);
end

% draw electrodes
if options.iselectrode
    electrcolor = [0.0  0.0  1.0];
    plot3(model.electrodes.locations(:,1), ...
        model.electrodes.locations(:,2), ...
        model.electrodes.locations(:,3), ...
        'k.','LineWidth',4,'color', electrcolor);
end

% draw labels
if options.islabel
    labelscale = 0.5;
    if isequal(options.channels,'single')
        if isequal(options.dispmodes, 'out2in')
            rr = options.whichchannel(1);
            cc = options.whichchannel(2);
            
            % the first location
            location = labelscale*model.electrodes.locations(rr,:);
            text( location(1), location(2), location(3), ...
                char(upper(model.electrodes.labels{rr})),'FontSize',8 ,...
                'HorizontalAlignment','center');
            
            % the second location
            location = labelscale*model.electrodes.locations(cc,:);
            text( location(1), location(2), location(3), ...
                char(upper(model.electrodes.labels{cc})),'FontSize',8 ,...
                'HorizontalAlignment','center');
        else
            i = options.whichchannel;
            location = labelscale*model.electrodes.locations(i,:);
            text( location(1), location(2), location(3), ...
                char(upper(model.electrodes.labels{i})),'FontSize',8 ,...
                'HorizontalAlignment','center');
        end
    else
        electrnum = length(model.electrodes.labels);
        for i = 1:electrnum
            location = labelscale*model.electrodes.locations(i,:);
            text( location(1), location(2), location(3), ...
                char(upper(model.electrodes.labels{i})),'FontSize',8 ,...
                'HorizontalAlignment','center');
        end
    end
end

% axis manual;
hfig = gcf;
if options.dtf
    roipos = model.electrodes.locations;
    
    ArrowSizeLimit = [1 10]*(options.displimits(2));
    SphereSizeLimit = [1 30];
    opt = struct('Channels', options.channels,...
        'Whichchannel', options.whichchannel,...
        'ValLim', options.displimits,...
        'ArSzLt',ArrowSizeLimit,...
        'SpSzLt',SphereSizeLimit);
    
    if isequal(options.dispmodes, 'out2in')
        drawdtfconngraph(currentdtfmatrix,roipos,opt);
    elseif isequal(options.dispmodes, 'outflow')
        outflows = sum(currentdtfmatrix,1)/(size(currentdtfmatrix,2)-1);
        drawdtfflowgraph(outflows,roipos,opt);
    else % inflow
        inflows = sum(currentdtfmatrix,2)/(size(currentdtfmatrix,2)-1);
        drawdtfflowgraph(inflows,roipos,opt);
    end
end

lighting phong; % phone, gouraud
lightcolor = [0.3 0.3 0.3];
light('Position',[0 0 1],'color',lightcolor);
light('Position',[0 1 0],'color',lightcolor);
light('Position',[0 -1 0],'color',lightcolor);
axis off;
