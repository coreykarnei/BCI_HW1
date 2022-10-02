clear
clc
load("selectedChannels.mat")
addpath(genpath('./toolboxes/'));
%loaded data from session 1
%[s,h]=sload('Subject_003*.gdf');
[s,h]=sload('Subject_003_TESS_Online__feedback__s001_r001_2021_07_06_153655.gdf');
%Only first 32 channels contain EEG data
s_EEG = s(:,1:32);


%alpha and beta bands
f_a = [9 11];
f_b = [18 22]; 
%number of samples
fs = 512;
%filtering to get just the mu band
N = 2;
[B, A] = butter(N, [f_a(1) f_b(2)]*2 / fs);
s_EEG_mu = filter(B, A, s_EEG);


total_trial_count = 1;
%looping through all events to extract important info
for i=1:numel(h.EVENT.TYP)
    rh = h.EVENT.TYP == 770; %for rh
    lh = h.EVENT.TYP == 769; %for lh
 
    trial_start = h.EVENT.TYP == 1000;
    trial_end1 = or(h.EVENT.TYP==7692,h.EVENT.TYP==7702);
    trial_end2 = or(h.EVENT.TYP==7693,h.EVENT.TYP==7703);
    trial_end = or(trial_end1,trial_end2);


    %creating table of times of event starts and ends
        if(trial_start(i)==1)
            trial_times(total_trial_count,1) = h.EVENT.POS(i);
        end
        if(trial_end(i)==1)
            trial_times(total_trial_count,2) = h.EVENT.POS(i);
            trial_times(total_trial_count,3) = h.EVENT.TYP(i);
            total_trial_count = total_trial_count + 1;
        end
end

total_trials = sum(trial_start);
disp(["total trials =",total_trials])



RH = sum(rh);
LH = sum(lh);

%before CAR filter
s_2 = s_EEG_mu.^2;

for j=1:size(s_EEG_mu,2)
    %doing CAR spatial filtering
    s_EEG_mu_car(:,j) = s_EEG_mu(:,j) - mean(s_EEG_mu,2);
end


%after CAR filter
s_2_mu = s_EEG_mu_car.^2;


%plotting s^2 for a specific channel(column)
%for miss using time 14111

    %file_title = '';
    cur_trial_num = 1;
    trial_num_7692 = 0;
    trial_num_7702 = 0;
    trial_num_7693 = 0;
    trial_num_7703 = 0;



for k=1:20
    t = tiledlayout(3,1);
    
    end_of_trial = trial_times(k,2);

    %Plotting s^2 all 32 channels in mu band
    nexttile
    plot(s_2_mu(end_of_trial-256:end_of_trial,:));
     title('Plot of all 32 channels')

    %Plotting grand average s^2 in mu band
    nexttile
    plot(mean(s_2_mu(end_of_trial-256:end_of_trial,:),2));
    title('Plot of grand average')
    
    %Topoplot of grand average s^2 in mu band
    nexttile
    topoplot_(mean(s_2_mu(end_of_trial-256:end_of_trial,:),2),selectedChannels);
    %colorbar('Ticks',min(color_ticks):(max(color_ticks)-min(color_ticks))/10:max(color_ticks))
    title('mu power')
    
    %creating files with appropriate names/trial numbers
    file_title = '';
%     cur_trial_num = 1;
%     trial_num_7692 = 0;
%     trial_num_7702 = 0;
%     trial_num_7693 = 0;
%     trial_num_7703 = 0;

    if(trial_times(k,3)==7692)
        trial_num_7692 = trial_num_7692+ 1;
        file_title = ['Trial ',num2str(cur_trial_num),' LH MISS_TIMEOUT ',num2str(trial_num_7692),'.pdf'];
        cur_trial_num = cur_trial_num + 1;
        exportgraphics(t,file_title);
    end
    
    if(trial_times(k,3)==7702)
        trial_num_7702 = trial_num_7702 + 1;
        file_title = ['Trial ',num2str(cur_trial_num),' RH MISS_TIMEOUT ',num2str(trial_num_7702),'.pdf'];
        cur_trial_num = cur_trial_num + 1;
        exportgraphics(t,file_title);
    end

    if(trial_times(k,3)==7693)
        trial_num_7693 = trial_num_7693+ 1;
        file_title = ['Trial ',num2str(cur_trial_num),' LH HIT ',num2str(trial_num_7693),'.pdf'];
        cur_trial_num = cur_trial_num + 1;
        exportgraphics(t,file_title);
    end

    if(trial_times(k,3)==7703)
        trial_num_7703 = trial_num_7703 + 1;
        file_title = ['Trial ',num2str(cur_trial_num),' RH HIT ',num2str(trial_num_7703),'.pdf'];
        cur_trial_num = cur_trial_num + 1;
        exportgraphics(t,file_title);
    end   
end














