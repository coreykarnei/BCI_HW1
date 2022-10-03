function [position] = subplotPosition(n,m,left,right,down,top,offsetx,offsety,ratiox,ratioy)

for i = 1:n
    for j = 1:m
        position(i,j,1) = (j-1)*ratiox/m +1*ratiox/m*left + offsetx;
        position(i,j,2) = (n-i)*ratioy/n +1*ratioy/n*down + offsety;
%         disp(['x = ' num2str((j-1)/m +1/m*xlabel) ',y = ' num2str((n-i)/n +1/n*ylabel)]);
        position(i,j,3) = 1*ratiox/m*(1-left-right);
        position(i,j,4) = 1*ratioy/n*(1-down-top);
    end
end