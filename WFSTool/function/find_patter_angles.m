function [alfa] = find_patter_angles(vsx,vsy,dir)
[x,y] = find_middle();
alfa = zeros(1,size(vsx,2));
for i = 1:size(vsx,2)
    alfa(i) = atan2(x - vsx(i),y - vsy(i));
    alfa(i) = -(alfa(i) .* 180 / pi +90);
    alfa(i) = alfa(i) - dir;
    if alfa(i) <= -180
        alfa(i) = alfa(i) + 360;
    end
end

