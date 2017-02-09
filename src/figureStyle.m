classdef figureStyle < handle
   
    properties 
        FontName = 'Arial';
        FontSize = 8
        DefaultColors = [
            9  63 112;        % Blue (0.21)
           16 153   9;        % Green (0.38)
          255 189 168]/255;   % Red (0.81)
        DefaultLineType = {'-','--',':'};
        DefaultMarkers = {'o','s','^'};
        DefaultLineLineWidth = 1.5;
        DefaultLineMarkerSize = 8;
    end
    
    properties (SetAccess = private)
        figureSize_names   = cell(0,1);
        figureSize_widths  = zeros(0,1);
        figureSize_heights = zeros(0,1);
    end
    
    methods
        function obj = figureStyle(style)
            switch lower(style)
                case 'asce'
                    
                case {'main','thesis'}
                    obj.addFigureSize('1x1',6.5,8.8);
                    obj.addFigureSize('2x1',6.5,4.3);
                    obj.addFigureSize('2x1s',4.3,4.3);
                    obj.addFigureSize('3x1',6.5,2.9);
                    obj.addFigureSize('3x2',3.1,2.9);
                
                case 'display'
                    obj.FontSize = 12;
                    
                otherwise
                    error('Unknown style: %s',style)
            end          
        end
        function addFigureSize(obj,name,width,height)
            assert(ischar(name),'Name should be a character string')
            assert(~ismember(name,obj.figureSize_names),'Figure size with name: %s already exists',name)
            obj.figureSize_names = vertcat(obj.figureSize_names,name);
            assert(isnumeric(width) && isscalar(width),'width should be a numeric scalar')
            obj.figureSize_widths = vertcat(obj.figureSize_widths,width);
            assert(isnumeric(height) && isscalar(height),'height should be a numeric scalar')
            obj.figureSize_heights = vertcat(obj.figureSize_heights,height);
        end
        function [width,height] = getFigureSize(obj,name)
            ind = find(strcmp(name,obj.figureSize_names));
            assert(~isempty(ind),'Unknown figure size name: %s',name)
            width  = obj.figureSize_widths(ind);
            height = obj.figureSize_heights(ind);
        end
        
        
        function h = figure(obj,width,height)
            % obj.figure(figureSize)
            % obj.figure(width,height)
            
            % Create figure
            h = figure;          
            set(h,'Color',[1 1 1]);
            set(h,'InvertHardcopy','off');
            
            if nargin > 1
                if ischar(width)
                    [width,height] = obj.getFigureSize(width);
                end
                set(h,'PaperUnits'       ,'inches');
                set(h,'PaperPositionMode','manual');
                set(h,'PaperSize'        ,[width height]);
                set(h,'PaperPosition'    ,[0 0 width height]);
                ppi = get(0,'ScreenPixelsPerInch');
                set(h,'Position',[100 100 width*ppi height*ppi]);
            end
            
            if nargout < 1
                clear h;
            end
        end           
        function h = axes(obj,position)
            % Create axis
            h = axes;

            % Set position
            if nargin > 1
                set(h,'Position',position);
            end
           
            % Defaults
            
            % Axes Fond
            set(h,'FontName',obj.FontName);
            set(h,'FontSize',obj.FontSize);
            
            % X and Y label font
            l = get(h,'XLabel');
            set(l,'FontName',obj.FontName);
            set(l,'FontSize',obj.FontSize);
            l = get(h,'YLabel');
            set(l,'FontName',obj.FontName);
            set(l,'FontSize',obj.FontSize);
            
            % 
            set(h,'DefaultLineLineWidth',obj.DefaultLineLineWidth);
            set(h,'DefaultLineMarkerSize',obj.DefaultLineMarkerSize);
            set(h,'TickLength',[0.02 0.02]);
            set(h,'Box','on');
            hold all
            
            if nargout < 1
                clear h;
            end
        end
        function h = backgroundAxes(obj)
            h = obj.axes([0 0 1 1]);
            set(gca,'visible','off');
        end
        function [hf,ha] = picture(obj,varargin)
            hf = obj.figure(varargin{:});
            ha = obj.backgroundAxes;
            axis equal
            if nargout < 2
                clear ha;
            end
            if nargout < 1
                clear hf;
            end
        end
    end
end




function setPaperSize(h,width,height,units)

end
