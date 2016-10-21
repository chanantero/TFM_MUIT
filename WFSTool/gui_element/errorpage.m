function [ ] = errorpage(str)
% The str here shows the error information
    fprintf(str);
    h = dialog('name','Error! ','position',[640,480,200,70]);  
    uicontrol('parent',h,'style','text','string',str,'position',[15 10 180 50],'fontsize',10);  
    uicontrol('parent',h,'style','pushbutton','string','OK','position',[75 5 50 30],'callback','delete(gcbf)');
end

