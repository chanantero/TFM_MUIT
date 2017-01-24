classdef WFSTool2 < handle
    
    properties
        fig
        player
    end
    
    methods
        
        function obj = WFSTool2()
            obj.player = reproductor();
            obj.fig = obj.createFigure();
%             obj.fig = openfig('WFSTool.fig');
%             
%             % Modify original callbacks
%             butTest = findobj(obj.fig.Children, 'Tag', 'but_test');
%             butTest.Callback = @(hObject, eventdata) obj.callback(); %@(hObject, eventdata) but_test_Callback(obj.playObj, Fs);
%             
%             butImportWav = findobj(obj.fig.Children, 'Tag', 'file');
%             table = findobj(obj.fig.Children, 'Tag', 'table');
%             butImportWav.Callback = @(hObject, eventdata) open_Callback(table);
        end
        
        function fig = createFigure(obj)
            fig = figure;
            
            reprodPanel = reproductionPanel(fig);
            
            
        end
        
        
        function playMusic(obj)
            % Audio file information
            fileName = 'Dangerous Woman.mp3';
            obj.player.audioFileName = fileName;
           
            % Reproduce
            obj.player.executeOrder('play');        
        end
       
        
        
    end
    
end

