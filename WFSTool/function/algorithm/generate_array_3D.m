function [coord, orient] = generate_array_3D()
% Generates the position and orientation of speakers in 2D.
% Output arguments:
% - coord. Nx3 array where N is the number of loudspeakers. Each of the
% three rows correspond to the x, y and z coordinates respectively.
% - orient. Nx3. Orientation of the loudspeakers. Each row contains the
% unitary orientation vector of the loudspeaker (main direction of the
% loudspeaker, broadside).
coord(1,1:8)=1.0182-(0.09:0.18:1.35)*cos(45*pi/180);
coord(2,1:8)=(0.09:0.18:1.35)*cos(45*pi/180);
coord(1,9:32)=0;
coord(2,9:32)=(1.0182)+(0.09:0.18:2*1.44+1.35);
coord(1,33:40)=(0.09:0.18:1.35)*cos(45*pi/180);
coord(2,33:40)=(1.0182+3*1.44)+(0.09:0.18:1.35)*cos(45*pi/180);
coord(1,41:48)=(1.0182)+(0.09:0.18:1.35);
coord(2,41:48)=(2*1.0182+3*1.44);
coord(1,49:56)=(1.0182+1.44)+(0.09:0.18:1.35)*cos(45*pi/180);
coord(2,49:56)=(2*1.0182+3*1.44)-(0.09:0.18:1.35)*cos(45*pi/180);
coord(1,57:80)=(2*1.0182+1.44);
coord(2,57:80)=(1.0182)+(2*1.44+1.35:-0.18:0.09);
coord(1,81:88)=(2*1.0182+1.44)-(0.09:0.18:1.35)*cos(45*pi/180);
coord(2,81:88)=(1.0182)-(0.09:0.18:1.35)*cos(45*pi/180);
coord(1,89:96)=(1.0182+1.44)-(0.09:0.18:1.35);
coord(2,89:96)=0;

% From 2D to 3D
N = size(coord', 1);
coord = [coord', zeros(N, 1)];

theta = [ones(1,8)*135,ones(1,24)*90,ones(1,8)*45,ones(1,8)*0,ones(1,8)*-45,ones(1,24)*-90,ones(1,8)*-135,ones(1,8)*180] - 90;
theta = theta';
orient = [cosd(theta), sind(theta), zeros(N, 1)];

end