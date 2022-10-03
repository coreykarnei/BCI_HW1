function eeg1=EEGFilterOrder(eeg,band,samplerate,order)

for i=1:size(eeg,1)
    if band(1)==0
        Wn=band(2)*2/samplerate;
        [a,b]=butter(order,Wn,'low');
    else
        Wn=[band(1) band(2)]*2/samplerate;
        [a,b]=butter(order,Wn,'bandpass');
    end
    eeg1(i,:)=filter(a,b,eeg(i,:));
end
