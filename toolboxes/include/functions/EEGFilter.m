function eeg1=EEGFilter(eeg,band,samplerate)

for i=1:size(eeg,1)
    if band(1)==0
        Wn=band(2)*2/samplerate;
        [a,b]=butter(2,Wn,'low');
    else
        Wn=[band(1) band(2)]*2/samplerate;
        [a,b]=butter(2,Wn,'bandpass');
    end
    eeg1(i,:)=filter(a,b,eeg(i,:));
end