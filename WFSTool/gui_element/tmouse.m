function tmouse(action)
global x y plotTimer;
 switch(action)
% Actions when first pressed
  case 'down',
   currPt = get(gca, 'CurrentPoint');
   if currPt(1,1) >=-2 && currPt(1,1) <=6 && currPt(1,2) >= -2 && currPt(1,2) <= 8
       x = currPt(1,1);
       y = currPt(1,2);
       set(gcf, 'WindowButtonMotionFcn', 'tmouse move');
       set(gcf, 'WindowButtonUpFcn', 'tmouse up');
       start(plotTimer);
   end

 % Actions when moving the mouse
  case 'move',
   currPt = get(gca, 'CurrentPoint');
   if currPt(1,1) >=-2 && currPt(1,1) <=6 && currPt(1,2) >= -2 && currPt(1,2) <= 8
       x = currPt(1,1);
       y = currPt(1,2);
   end
% Actions when release the mouse
  case 'up',
   set(gcf, 'WindowButtonMotionFcn', '');
   set(gcf, 'WindowButtonUpFcn', '');
   stop(plotTimer);
 end