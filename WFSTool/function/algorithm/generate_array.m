function [alt] = generate_array()
% Output arguments:
% - alt. 2xN array where N is the number of loudspeakers. The first
% component of each column contains the x coordinate of the loudspeaker and
% the second component contains the y coordinate.
alt(1,1:8)=1.0182-(0.09:0.18:1.35)*cos(45*pi/180);
alt(2,1:8)=(0.09:0.18:1.35)*cos(45*pi/180);
alt(1,9:32)=0;
alt(2,9:32)=(1.0182)+(0.09:0.18:2*1.44+1.35);
alt(1,33:40)=(0.09:0.18:1.35)*cos(45*pi/180);
alt(2,33:40)=(1.0182+3*1.44)+(0.09:0.18:1.35)*cos(45*pi/180);
alt(1,41:48)=(1.0182)+(0.09:0.18:1.35);
alt(2,41:48)=(2*1.0182+3*1.44);
alt(1,49:56)=(1.0182+1.44)+(0.09:0.18:1.35)*cos(45*pi/180);
alt(2,49:56)=(2*1.0182+3*1.44)-(0.09:0.18:1.35)*cos(45*pi/180);
alt(1,57:80)=(2*1.0182+1.44);
alt(2,57:80)=(1.0182)+(2*1.44+1.35:-0.18:0.09);
alt(1,81:88)=(2*1.0182+1.44)-(0.09:0.18:1.35)*cos(45*pi/180);
alt(2,81:88)=(1.0182)-(0.09:0.18:1.35)*cos(45*pi/180);
alt(1,89:96)=(1.0182+1.44)-(0.09:0.18:1.35);
alt(2,89:96)=0;
end

