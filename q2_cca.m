load('ErrP_data_scripts\ErrP_channels.mat')
load('ErrP_data_scripts\ErrP_data_HW1.mat')
c0=0;
Y1=zeros(32,1024);
Y2=zeros(32,1024);
yy = permute(trainingEpochs.rotation_data,[2 1 3]);
for i=1:320
    if trainingEpochs.label(i) == 0
        Y1 = cat(3,Y1,yy(:,:,i));
        c0 = c0+1;
    else
        Y2 = cat(3,Y2,yy(:,:,i));
    end
end
Y1 = Y1(:,:,2:end);
Y2 = Y2(:,:,2:end);
Y = cat(3,Y1,Y2);

Y = reshape(Y, [32, 1024*320]);

X1 = repmat(mean(Y1,3),1,c0);
X2 = repmat(mean(Y2,3),1,320-c0);
X = [X1 X2];

[A,B,R,U,V] = canoncorr(X,Y);

save('u.mat','U')