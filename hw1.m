[s,h]=sload('Subject_003*.gdf');
%alpha and beta bands
f_a = [9 11];
f_b = [18 22]; 
f_mu = [8,12];
%number of samples
fs = 512;
%filtering to get just the mu band
N = 2;
[B, A] = butter(N, [f_mu(1) f_mu(2)]*2 / fs);
s_mu = filter(B, A, s);
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
            trial_times(i,1) = h.EVENT.POS(i);
        end
    
        if(trial_end(i)==1)
            trial_times(i,1) = h.EVENT.POS(i);
        end
end
temp = 1;
for l=1:numel(trial_times)
    if(trial_times(l)~=0)
        if(mod(temp,2)~=0)
            final_times(temp,1) = trial_times(l);
            t_start_time = trial_times(l);
        else
            final_times(temp-1,2) = trial_times(l);
            t_end_time = trial_times(l);
        end
    temp = temp +1;    
    end 
end
% for z=1:numel(final_times)
%     if(final_times(z)~=0)
%        trial_s2 = s_mu_reduced(final_times(z,1):final_times(z,2),:); 
%     end
% end
total_trials = sum(trial_start);
RH = sum(rh);
LH = sum(lh);
s_mu_reduced = s_mu(:,1:32);
for j=1:size(s_mu_reduced,2)
    %doing CAR spatial filtering
    s_mu_reduced(:,j) = s_mu_reduced(:,j) - mean(s_mu_reduced,2);
end
s_2 = s_mu_reduced.^2;
%l = length(s_2(14111-256,:):s_2(14111,:)));
plot(s_2(14111-256:14111,2));

disp(s(1,1));
s11_squared = s(1,1).^2;
disp(s11_squared)
disp(s_2(1,1));