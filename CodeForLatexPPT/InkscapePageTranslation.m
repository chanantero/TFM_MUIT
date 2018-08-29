refPos1 = [0 350.634];
refPos2 = [18.504 -812.621];

pagePos = [0 0 980.407 1074.226]; % x y width height

desp = refPos2 - refPos1;

newPagePos = pagePos;
newPagePos(1:2) = newPagePos(1:2) + desp;

refPos1 - pagePos(1:2)
refPos2 - newPagePos(1:2)
