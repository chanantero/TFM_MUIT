function [ ] = PlayTimer(obj,event)
% The timer called for playing the music.
global music pageSize startSample endSample endPoint pageNumList;
global parray_act an tn isActive playTimer amp_z pattern angles;
fileCount = size(music,2);
page_no = fileCount * 2;
pageBufCount = 10;
runMaxSpeed = false;
oldlist = [];
totallist = [];
max_endPoint = max(endPoint);
amp = 1;

for a = 1:page_no
    oldpage = [];
    % If the music is all finished
    if startSample > max_endPoint
        stop(playTimer);
        fprintf('Music Finished. \n');
        return;
    end
    for b = 1:fileCount
        ts = startSample;
        if ts < endPoint(b) && isActive(b) == 1
            % For each source, calculate pattern attenuation first
            for o = 1:size(pattern,2)
                if pattern(1,o) == round(angles(b))
                    amp = pattern(2,o);
                    break;
                end
            end
            att = amp * amp_z;
            % If the attenuation is zero, jump the source
            if att ~= 0
                te = endSample; 
                k = cell2mat(music(b));
                x = k(ts-pageSize:te-pageSize);
                if te > endPoint(b)
                    te = endPoint(b);
                end
                y = k(ts:te);
                y = [y;zeros(pageSize-size(y,1),1)];
                newlist = cell2mat(parray_act(b));
                newpage = WFS_implement(y,newlist,an(b,:),tn(b,:),pageSize,x);
                newpage = att * newpage;
                if isempty(oldpage)
                    totallist = newlist;
                    oldlist = newlist;
                    oldpage = newpage;
                else
                    totallist = [oldlist,newlist];
                    Or_old = zeros(pageSize,96);
                    Or_new = zeros(pageSize,96);
                    for i = 1:size(oldlist,2)
                        c = oldlist(i);
                        Or_old(:,c) = oldpage(:,i);
                    end
                    for i = 1:size(newlist,2)
                        c = newlist(i);
                        Or_new(:,c) = newpage(:,i);
                    end   
                    oldpage = Or_old + Or_new;
                    oldpage(:,find(sum(abs(oldpage),1) == 0)) = [];
                    totallist = unique(totallist);
                    oldlist = totallist;
                end
            end
        end  
    end
    
    if isempty(totallist) == 1
       stop(playTimer);
       errorpage('No valid sound source! ');
       return;
    end
    
    if size(oldpage,2) < size(totallist,2)
        oldpage = zeros(pageSize,size(totallist,2));
    end
    
    pageNumList = [pageNumList playrec('play',oldpage,totallist)];
    if(startSample == 1)
        playrec('resetSkippedSampleCount');
    end
    if(length(pageNumList) > pageBufCount)
        if(runMaxSpeed)
            while(playrec('isFinished', pageNumList(1)) == 0)
            end
        else        
            playrec('block', pageNumList(1));
        end
        pageNumList = pageNumList(2:end);
    end
    startSample = startSample + pageSize;
    endSample = endSample + pageSize;  
end
end

