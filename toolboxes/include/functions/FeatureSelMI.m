% Feature selection with mutual information
function index = FeatureSelMI(label,feature,num)

for i = 1:size(feature,1)
    feature(i,:) = feature(i,:)/max(feature(i,:));
end

index = [];
for i = 1:num
    hist(feature(1,:));
    for j = 1:size(feature,1)
        entropy(feature(j,:));
%         if i == 1
%             I(j) = entropy(label) + entropy(feature(j,:)) - jentropy([label;feature(j,:)]);
%         else
            I(j) = entropy(label) + jentropy(feature([index j],:)) - jentropy([label;feature([index j],:)]);
%         end
    end
    [mvalue mindex] = max(I);
    index = [index mindex];
end