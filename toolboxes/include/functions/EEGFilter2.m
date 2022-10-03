function eeg2 = EEGFilter2(eeg,band,samplerate)

for i=1:size(eeg,1)
    if band(1)==0
        Wn=band(2)*2/samplerate;
        [a,b]=butter(4,Wn,'low');
        eeg2(i,:)=filter(a,b,eeg(i,:));
    else
        Wn=band(2)*2/samplerate;
        [a1,b1]=butter(4,Wn,'low');
        
        Wn=band(1)*2/samplerate;
        [a2,b2]=butter(4,Wn,'high');
        
        eeg1(i,:)=filter(a1,b1,eeg(i,:));
        eeg2(i,:)=filter(a2,b2,eeg1(i,:));
    end
end