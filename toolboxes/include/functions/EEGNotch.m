function eeg = EEGNotch(eeg,freq,fs);

flag = 0;
if size(eeg,1) > size(eeg,2)
    eeg = eeg';
    flag = 1;
end

for i=1:size(eeg,1)
    for j=1:length(freq)
        wo = freq(j)/(fs/2);  bw = wo/20;
        [b,a] = iirnotch(wo,bw);
        eeg(i,:)=filter(b,a,eeg(i,:));
    end
end

if flag == 1
    eeg = eeg';
end