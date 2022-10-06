clear
clc
load("selectedChannels.mat")
addpath(genpath('./toolboxes/'));
%loaded data from sessions 1-4
[s,h]=sload('Subject_003_TESS_Online__feedback__s001_r001_2021_07_06_153655.gdf');
[s2,h2]=sload('Subject_003_TESS_Online__feedback__s001_r002_2021_07_06_154358.gdf');
[s3,h3]=sload('Subject_003_TESS_Online__feedback__s001_r003_2021_07_06_155155.gdf');
[s4,h4]=sload('Subject_003_TESS_Online__feedback__s001_r004_2021_07_06_160446.gdf');
%Only first 32 channels contain EEG data
s_EEG = s(:,1:32);
s_EEG2 = s2(:,1:32);
s_EEG3 = s3(:,1:32);
s_EEG4 = s4(:,1:32);


%alpha and beta bands
%f_a = [9 13];
f_a = [8 12];
f_b = [18 22]; 
%number of samples
fs = 512;
%filtering to get just the mu band
N = 2;
[B, A] = butter(N, [f_a(1) f_a(2)]*2 / fs,"bandpass");
s_EEG_mu = filtfilt(B, A, s_EEG);
s_EEG_mu2 = filtfilt(B, A, s_EEG2);
s_EEG_mu3 = filtfilt(B, A, s_EEG3);
s_EEG_mu4 = filtfilt(B, A, s_EEG4);


total_trial_count = 1;
%looping through all events to extract important info
for i=1:numel(h.EVENT.TYP)
    trial_times(total_trial_count,4) = 1;
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

for i=1:numel(h2.EVENT.TYP)
    trial_times(total_trial_count,4) = 2;
    rh = h2.EVENT.TYP == 770; %for rh
    lh = h2.EVENT.TYP == 769; %for lh
 
    trial_start = h2.EVENT.TYP == 1000;
    trial_end1 = or(h2.EVENT.TYP==7692,h2.EVENT.TYP==7702);
    trial_end2 = or(h2.EVENT.TYP==7693,h2.EVENT.TYP==7703);
    trial_end = or(trial_end1,trial_end2);
    %creating table of times of event starts and ends
        if(trial_start(i)==1)
            trial_times(total_trial_count,1) = h2.EVENT.POS(i);
        end
        if(trial_end(i)==1)
            trial_times(total_trial_count,2) = h2.EVENT.POS(i);
            trial_times(total_trial_count,3) = h2.EVENT.TYP(i);
            total_trial_count = total_trial_count + 1;
        end
end

for i=1:numel(h3.EVENT.TYP)
    trial_times(total_trial_count,4) = 3;
    rh = h3.EVENT.TYP == 770; %for rh
    lh = h3.EVENT.TYP == 769; %for lh
 
    trial_start = h3.EVENT.TYP == 1000;
    trial_end1 = or(h3.EVENT.TYP==7692,h3.EVENT.TYP==7702);
    trial_end2 = or(h3.EVENT.TYP==7693,h3.EVENT.TYP==7703);
    trial_end = or(trial_end1,trial_end2);
    %creating table of times of event starts and ends
        if(trial_start(i)==1)
            trial_times(total_trial_count,1) = h3.EVENT.POS(i);
        end
        if(trial_end(i)==1)
            trial_times(total_trial_count,2) = h3.EVENT.POS(i);
            trial_times(total_trial_count,3) = h3.EVENT.TYP(i);
            total_trial_count = total_trial_count + 1;
        end
end

for i=1:numel(h4.EVENT.TYP)
    trial_times(total_trial_count,4) = 4;
    rh = h4.EVENT.TYP == 770; %for rh
    lh = h4.EVENT.TYP == 769; %for lh
 
    trial_start = h4.EVENT.TYP == 1000;
    trial_end1 = or(h4.EVENT.TYP==7692,h4.EVENT.TYP==7702);
    trial_end2 = or(h4.EVENT.TYP==7693,h4.EVENT.TYP==7703);
    trial_end = or(trial_end1,trial_end2);
    %creating table of times of event starts and ends
        if(trial_start(i)==1)
            trial_times(total_trial_count,1) = h4.EVENT.POS(i);
        end
        if(trial_end(i)==1)
            trial_times(total_trial_count,2) = h4.EVENT.POS(i);
            trial_times(total_trial_count,3) = h4.EVENT.TYP(i);
            total_trial_count = total_trial_count + 1;
        end
end

total_trials = size(trial_times)-1;
disp(["total trials =",total_trials])
LHMissTimeout=1;
RHMissTimeout=1;
LHHit=1;
RHHit=1;
Total_LH = 1;
for trial=1:80
    %Class LH MISS/TIMEOUT
    if(trial_times(trial,3)==7692)
        LHMissTimeout_table(LHMissTimeout,1) = trial_times(trial,1);
        LHMissTimeout_table(LHMissTimeout,2) = trial_times(trial,2);
        LHMissTimeout_table(LHMissTimeout,3) = trial_times(trial,3);

        if(trial_times(trial,4)==1) 
            LHMissTimeout_table(LHMissTimeout,4)=1;
            Total_LH_t(Total_LH,4)=1;
        end

        if(trial_times(trial,4)==2) 
            LHMissTimeout_table(LHMissTimeout,4)=2;
            Total_LH_t(Total_LH,4)=2;
        end

        if(trial_times(trial,4)==3) 
            LHMissTimeout_table(LHMissTimeout,4)=3;
            Total_LH_t(Total_LH,4)=3;
        end
        if(trial_times(trial,4)==4) 
            LHMissTimeout_table(LHMissTimeout,4)=4; 
            Total_LH_t(Total_LH,4)=4;        
        end
        LHMissTimeout = LHMissTimeout + 1;
        Total_LH = Total_LH + 1;
    end
    %Class RH MISS/TIMEOUT
    if(trial_times(trial,3)==7702)
        RHMissTimeout_table(RHMissTimeout,1) = trial_times(trial,1);
        RHMissTimeout_table(RHMissTimeout,2) = trial_times(trial,2);
        RHMissTimeout_table(RHMissTimeout,3) = trial_times(trial,3);
        if(trial_times(trial,4)==1) RHMissTimeout_table(RHMissTimeout,4)=1; end
        if(trial_times(trial,4)==2) RHMissTimeout_table(RHMissTimeout,4)=2; end
        if(trial_times(trial,4)==3) RHMissTimeout_table(RHMissTimeout,4)=3; end
        if(trial_times(trial,4)==4) RHMissTimeout_table(RHMissTimeout,4)=4; end
        RHMissTimeout = RHMissTimeout + 1;
    end
    %Class LH HIT
    if(trial_times(trial,3)==7693)
        LHHit_table(LHHit,1) = trial_times(trial,1);
        LHHit_table(LHHit,2) = trial_times(trial,2);
        LHHit_table(LHHit,3) = trial_times(trial,3);
        if(trial_times(trial,4)==1) LHHit_table(LHHit,4)=1; end
        if(trial_times(trial,4)==2) LHHit_table(LHHit,4)=2; end
        if(trial_times(trial,4)==3) LHHit_table(LHHit,4)=3; end
        if(trial_times(trial,4)==4) LHHit_table(LHHit,4)=4; end
        LHHit = LHHit + 1;
    end
    %Class RH HIT (apparently none?)
    if(trial_times(trial,3)==7703)
        RHHit_table(RHHit,1) = trial_times(trial,1);
        RHHit_table(RHHit,2) = trial_times(trial,2);
        RHHit_table(RHHit,3) = trial_times(trial,3); 
        RHHit = RHHit + 1;
    end

   


end

 LHTotal_table = [LHHit_table;LHMissTimeout_table];
 RHTotal_table = RHMissTimeout_table;

%before CAR filter
s_2 = s_EEG_mu.^2;
s_2_2 = s_EEG_mu2.^2;
s_2_3 = s_EEG_mu3.^2;
s_2_4 = s_EEG_mu4.^2;

for j=1:size(s_EEG_mu,2)
    %doing CAR spatial filtering
    s_EEG_mu_car(:,j) = s_EEG_mu(:,j) - mean(s_EEG_mu,2);
end
for j=1:size(s_EEG_mu,2)
    %doing CAR spatial filtering
    s_EEG_mu_car2(:,j) = s_EEG_mu2(:,j) - mean(s_EEG_mu2,2);
end

for j=1:size(s_EEG_mu,2)
    %doing CAR spatial filtering
    s_EEG_mu_car3(:,j) = s_EEG_mu3(:,j) - mean(s_EEG_mu3,2);
end

for j=1:size(s_EEG_mu,2)
    %doing CAR spatial filtering
    s_EEG_mu_car4(:,j) = s_EEG_mu4(:,j) - mean(s_EEG_mu4,2);
end



%after CAR filter
s_2_mu = s_EEG_mu_car.^2;
s_2_mu2 = s_EEG_mu_car2.^2;
s_2_mu3 = s_EEG_mu_car3.^2;
s_2_mu4 = s_EEG_mu_car4.^2;

for k=1:40
    end_of_trial = LHTotal_table(k,2);
    if(LHTotal_table(k,4)==1)
    LHTotal_mean_CAR(:,:,k) = s_2_mu(end_of_trial-256:end_of_trial,:);
    LHTotal_mean(:,:,k) = s_2(end_of_trial-256:end_of_trial,:);
    end
    if(LHTotal_table(k,4)==2)
    LHTotal_mean_CAR(:,:,k) = s_2_mu2(end_of_trial-256:end_of_trial,:);
    LHTotal_mean(:,:,k) = s_2_2(end_of_trial-256:end_of_trial,:);
    end
    if(LHTotal_table(k,4)==3)
    LHTotal_mean_CAR(:,:,k) = s_2_mu3(end_of_trial-256:end_of_trial,:);
    LHTotal_mean(:,:,k) = s_2_3(end_of_trial-256:end_of_trial,:);
    end
    if(LHTotal_table(k,4)==4)
    LHTotal_mean_CAR(:,:,k) = s_2_mu4(end_of_trial-256:end_of_trial,:);
    LHTotal_mean(:,:,k) = s_2_4(end_of_trial-256:end_of_trial,:);
    end
end

for k=1:40
    end_of_trial = RHTotal_table(k,2);
    if(RHTotal_table(k,4)==1)
    RHTotal_mean_CAR(:,:,k) = s_2_mu(end_of_trial-256:end_of_trial,:);
    RHTotal_mean(:,:,k) = s_2(end_of_trial-256:end_of_trial,:);
    end
    if(RHTotal_table(k,4)==2)
    RHTotal_mean_CAR(:,:,k) = s_2_mu2(end_of_trial-256:end_of_trial,:);
    RHTotal_mean(:,:,k) = s_2_2(end_of_trial-256:end_of_trial,:);
    end
    if(RHTotal_table(k,4)==3)
    RHTotal_mean_CAR(:,:,k) = s_2_mu3(end_of_trial-256:end_of_trial,:);
    RHTotal_mean(:,:,k) = s_2_3(end_of_trial-256:end_of_trial,:);
    end
    if(RHTotal_table(k,4)==4)
    RHTotal_mean_CAR(:,:,k) = s_2_mu4(end_of_trial-256:end_of_trial,:);
    RHTotal_mean(:,:,k) = s_2_4(end_of_trial-256:end_of_trial,:);
    end
end


for k=1:size(LHHit_table)
    end_of_trial = LHHit_table(k,2);
    if(LHHit_table(k,4)==1)
        LHHit_mean_CAR(:,k) = mean(s_2_mu(end_of_trial-256:end_of_trial,:),2);
        LHHit_mean(:,k) = mean(s_2(end_of_trial-256:end_of_trial,:),2);
    end
    if(LHHit_table(k,4)==2)
        LHHit_mean_CAR(:,k) = mean(s_2_mu2(end_of_trial-256:end_of_trial,:),2);
        LHHit_mean(:,k) = mean(s_2_2(end_of_trial-256:end_of_trial,:),2);
    end
    if(LHHit_table(k,4)==3)
        LHHit_mean_CAR(:,k) = mean(s_2_mu3(end_of_trial-256:end_of_trial,:),2);
        LHHit_mean(:,k) = mean(s_2_3(end_of_trial-256:end_of_trial,:),2);
    end
    if(LHHit_table(k,4)==4)
        LHHit_mean_CAR(:,k) = mean(s_2_mu4(end_of_trial-256:end_of_trial,:),2);
        LHHit_mean(:,k) = mean(s_2_4(end_of_trial-256:end_of_trial,:),2);
    end    
end

%for LHMissTimeout
for k=1:size(LHMissTimeout_table)
    end_of_trial = LHMissTimeout_table(k,2);
    if(LHMissTimeout_table(k,4)==1)
        LHMissTimeout_mean_CAR(:,k) = mean(s_2_mu(end_of_trial-256:end_of_trial,:),2);
        LHMissTimeout_mean(:,k) = mean(s_2(end_of_trial-256:end_of_trial,:),2);
    end
    if(LHMissTimeout_table(k,4)==2)
        LHMissTimeout_mean_CAR(:,k) = mean(s_2_mu2(end_of_trial-256:end_of_trial,:),2);
        LHMissTimeout_mean(:,k) = mean(s_2_2(end_of_trial-256:end_of_trial,:),2);
    end
    if(LHMissTimeout_table(k,4)==3)
        LHMissTimeout_mean_CAR(:,k) = mean(s_2_mu3(end_of_trial-256:end_of_trial,:),2);
        LHMissTimeout_mean(:,k) = mean(s_2_3(end_of_trial-256:end_of_trial,:),2);
    end
    if(LHMissTimeout_table(k,4)==4)
        LHMissTimeout_mean_CAR(:,k) = mean(s_2_mu4(end_of_trial-256:end_of_trial,:),2);
        LHMissTimeout_mean(:,k) = mean(s_2_4(end_of_trial-256:end_of_trial,:),2);
    end    
end

%RHMissTimeout
for k=1:size(RHMissTimeout_table)
    end_of_trial = RHMissTimeout_table(k,2);
    if(RHMissTimeout_table(k,4)==1)
        RHMissTimeout_mean_CAR(:,k) = mean(s_2_mu(end_of_trial-256:end_of_trial,:),2);
        RHMissTimeout_mean(:,k) = mean(s_2(end_of_trial-256:end_of_trial,:),2);
    end
    if(RHMissTimeout_table(k,4)==2)
        RHMissTimeout_mean_CAR(:,k) = mean(s_2_mu2(end_of_trial-256:end_of_trial,:),2);
        RHMissTimeout_mean(:,k) = mean(s_2_2(end_of_trial-256:end_of_trial,:),2);
    end
    if(RHMissTimeout_table(k,4)==3)
        RHMissTimeout_mean_CAR(:,k) = mean(s_2_mu3(end_of_trial-256:end_of_trial,:),2);
        RHMissTimeout_mean(:,k) = mean(s_2_3(end_of_trial-256:end_of_trial,:),2);
    end
    if(RHMissTimeout_table(k,4)==4)
        RHMissTimeout_mean_CAR(:,k) = mean(s_2_mu4(end_of_trial-256:end_of_trial,:),2);
        RHMissTimeout_mean(:,k) = mean(s_2_4(end_of_trial-256:end_of_trial,:),2);
    end    
end



color_ticks = [0 2];
color_ticksM = [6.5 10];

%LH Class
figure()
t = tiledlayout(2,1);
nexttile
plot(mean(LHTotal_mean_CAR, [2,3]));
title('Plot of LH grand average WITH CAR')

nexttile
LH_CAR_TOPO = mean(LHTotal_mean_CAR, [1,3]); 
topoplot_(LH_CAR_TOPO,selectedChannels,'electrodes','labels','maplimits', [min(color_ticks) max(color_ticks)]);
colorbar('Ticks',min(color_ticks):(max(color_ticks)-min(color_ticks))/10:max(color_ticks))
title('LH Grand Average of MU power WITH CAR')

figure()
t = tiledlayout(2,1);

nexttile
plot(mean(LHTotal_mean, [2,3]));
title('Plot of LH grand average WITHOUT CAR')

nexttile
LH_NO_CAR_TOPO = mean(LHTotal_mean, [1,3]); 
topoplot_(LH_NO_CAR_TOPO,selectedChannels,'electrodes','labels','maplimits', [min(color_ticksM) max(color_ticksM)]);
colorbar('Ticks',min(color_ticksM):(max(color_ticksM)-min(color_ticksM))/10:max(color_ticksM))
title('LH Grand Average of MU power WITHOUT CAR')

%RH Class
figure()
t = tiledlayout(2,1);
nexttile
plot(mean(RHTotal_mean_CAR, [2,3]));
title('Plot of RH grand average WITH CAR')

nexttile
RH_CAR_TOPO = mean(RHTotal_mean_CAR, [1,3]); 
topoplot_(RH_CAR_TOPO,selectedChannels,'electrodes','labels','maplimits', [min(color_ticks) max(color_ticks)]);
colorbar('Ticks',min(color_ticks):(max(color_ticks)-min(color_ticks))/10:max(color_ticks))
title('RH Grand Average of MU power WITH CAR')

figure()
t = tiledlayout(2,1);

nexttile
plot(mean(RHTotal_mean, [2,3]));
title('Plot of RH grand average WITHOUT CAR')

nexttile
RH_NO_CAR_TOPO = mean(RHTotal_mean, [1,3]); 
topoplot_(RH_NO_CAR_TOPO,selectedChannels,'electrodes','labels','maplimits', [min(color_ticksM) max(color_ticksM)]);
colorbar('Ticks',min(color_ticksM):(max(color_ticksM)-min(color_ticksM))/10:max(color_ticksM))
title('RH Grand Average of MU power WITHOUT CAR')


% 
% %LH Hit
% figure()
% t = tiledlayout(2,1);
% nexttile
% LHHit_mean_GA_CAR = mean(LHHit_mean_CAR,2);
% plot(LHHit_mean_GA_CAR);
% title('Plot of grand average (w/CAR)')
% 
% nexttile
% topoplot_(LHHit_mean_GA_CAR,selectedChannels,'electrodes','labels','maplimits', [min(color_ticks) max(color_ticks)]);
% colorbar('Ticks',min(color_ticks):(max(color_ticks)-min(color_ticks))/10:max(color_ticks))
% title('Grand Average of MU power for LH Hit (CAR)')
% 
% figure()
% t = tiledlayout(2,1);
% nexttile
% LHHit_mean_GA = mean(LHHit_mean,2);
% plot(LHHit_mean_GA);
% title('Plot of grand average (w/o CAR)')
% 
% nexttile
% topoplot_(LHHit_mean_GA,selectedChannels,'electrodes','labels','maplimits', [min(color_ticksM) max(color_ticksM)]);
% colorbar('Ticks',min(color_ticksM):(max(color_ticksM)-min(color_ticksM))/10:max(color_ticksM))
% title('Grand Average of MU power for LH Hit (no CAR)')
% 
% %LH MissTimeout
% figure()
% t = tiledlayout(2,1);
% nexttile
% LHMissTimeout_mean_GA_CAR = mean(LHMissTimeout_mean_CAR,2);
% plot(LHMissTimeout_mean_GA_CAR);
% title('Plot of grand average (w/CAR)')
% 
% nexttile
% topoplot_(LHMissTimeout_mean_GA_CAR,selectedChannels,'electrodes','labels','maplimits', [min(color_ticks) max(color_ticks)]);
% colorbar('Ticks',min(color_ticks):(max(color_ticks)-min(color_ticks))/10:max(color_ticks))
% title('Grand Average of MU power for LH Miss/Timeout (CAR)')
% 
% figure()
% t = tiledlayout(2,1);
% nexttile
% LHMissTimeout_mean_GA = mean(LHMissTimeout_mean,2);
% plot(LHMissTimeout_mean_GA);
% title('Plot of grand average (w/o CAR)')
% 
% nexttile
% topoplot_(LHMissTimeout_mean_GA,selectedChannels,'electrodes','labels','maplimits', [min(color_ticksM) max(color_ticksM)]);
% colorbar('Ticks',min(color_ticksM):(max(color_ticksM)-min(color_ticksM))/10:max(color_ticksM))
% title('Grand Average of MU power for LH Miss/Timeout (no CAR)')
% 
% %RH MissTimeout
% figure()
% t = tiledlayout(2,1);
% nexttile
% RHMissTimeout_mean_GA_CAR = mean(RHMissTimeout_mean_CAR,2);
% plot(RHMissTimeout_mean_GA_CAR);
% title('Plot of grand average (w/CAR)')
% 
% nexttile
% topoplot_(RHMissTimeout_mean_GA_CAR,selectedChannels,'electrodes','labels','maplimits', [min(color_ticks) max(color_ticks)]);
% colorbar('Ticks',min(color_ticks):(max(color_ticks)-min(color_ticks))/10:max(color_ticks))
% title('Grand Average of MU power for RH Miss/Timeout (CAR)')
% 
% figure()
% t = tiledlayout(2,1);
% nexttile
% RHMissTimeout_mean_GA = mean(RHMissTimeout_mean,2);
% plot(RHMissTimeout_mean_GA);
% title('Plot of grand average (w/o CAR)')
% 
% nexttile
% topoplot_(RHMissTimeout_mean_GA,selectedChannels,'electrodes','labels','maplimits', [min(color_ticksM) max(color_ticksM)]);
% colorbar('Ticks',min(color_ticksM):(max(color_ticksM)-min(color_ticksM))/10:max(color_ticksM))
% title('Grand Average of MU power for RH Miss/Timeout (no CAR)')






 
    














