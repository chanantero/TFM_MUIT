function [ image ] = image_rotate(angle,file)
% By inputting rotating angle and the name of the file, the function will 
% return the image utility. 
img = imread(file);
[h,w] = size(img);

theta = angle/180*pi;
rot=[cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1];
pix1=[1 1 1]*rot;
pix2=[1 w 1]*rot;
pix3=[h 1 1]*rot;
pix4=[h w 1]*rot;

height=round(max(abs(pix1(1)-pix4(1))+0.5,abs(pix2(1)-pix3(1))+0.5));
width=round(max(abs(pix1(2)-pix4(2))+0.5,abs(pix2(2)-pix3(2))+0.5));
imgn=zeros(height,width);

delta_x=abs(min([pix1(1) pix2(1) pix3(1) pix4(1)]));
delta_y=abs(min([pix1(2) pix2(2) pix3(2) pix4(2)]));

for i=1-delta_y:height-delta_y
    for j=1-delta_x:width-delta_x
        pix=[i j 1]*rot;
        float_y=pix(1)-floor(pix(1));
        float_x=pix(2)-floor(pix(2));
        if pix(1)>=1 && pix(2) >= 1 && pix(1) <= h && pix(2) <= w
            pix_up_left=[floor(pix(1)) floor(pix(2))];
            pix_up_right=[floor(pix(1)) ceil(pix(2))];
            pix_down_left=[ceil(pix(1)) floor(pix(2))];
            pix_down_right=[ceil(pix(1)) ceil(pix(2))];
            value_up_left=(1-float_x)*(1-float_y);
            value_up_right=float_x*(1-float_y);
            value_down_left=(1-float_x)*float_y;
            value_down_right=float_x*float_y;
            imgn(i+delta_y,j+delta_x)=value_up_left*img(pix_up_left(1),pix_up_left(2))+...
                value_up_right*img(pix_up_right(1),pix_up_right(2))+...
                value_down_left*img(pix_down_left(1),pix_down_left(2))+...
                value_down_right*img(pix_down_right(1),pix_down_right(2));
        end
    end
end
image = uint8(imgn);
end

