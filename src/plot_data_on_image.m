classdef plot_data_on_image < handle
    
    properties
        image_filename
        p1_cimg
        p1_cdata
        p2_cimg
        p2_cdata
        p3_cimg
        p3_cdata
    end
    
    methods
        %% Constructor
        function obj = plot_data_on_image(image_filename)
            obj.image_filename = image_filename;
        end
        
        %% Setting Functions
        function set_point(obj,point,image_coords,data_coords)
            switch point
                case 1
                    obj.p1_cimg  = image_coords;
                    obj.p1_cdata = data_coords;
                case 2
                    obj.p2_cimg  = image_coords;
                    obj.p2_cdata = data_coords;
                case 3
                    obj.p3_cimg  = image_coords;
                    obj.p3_cdata = data_coords;
                otherwise
                    error('Invalid point')
            end
        end
        
        %% Transformation Function
        function [xi,yi] = transform(obj,x,y)
            % Check input
            assert(isnumeric(x), 'x should be numeric')
            if isrow(x)
                x = x';
            end
            assert(iscolumn(x), 'x should be a column vector')
            assert(isnumeric(y), 'y should be numeric')
            if isrow(y)
                y = y';
            end
            assert(iscolumn(y), 'y should be a column vector')
            
            % Define Transformation
            A = [ 1 obj.p1_cdata(1) obj.p1_cdata(2) 
                  1 obj.p2_cdata(1) obj.p2_cdata(2)
                  1 obj.p3_cdata(1) obj.p3_cdata(2) ];
            cimgx = linsolve(A,[obj.p1_cimg(1) obj.p2_cimg(1) obj.p3_cimg(1)]');
            cimgy = linsolve(A,[obj.p1_cimg(2) obj.p2_cimg(2) obj.p3_cimg(2)]');
            
            % Transform and Plot Data
            xx = [ones(size(x)) x y];
            xi = xx*cimgx;
            yi = xx*cimgy;            
        end
        
        %% Plotting 
        function show_image(obj)
            imshow(obj.image_filename,'InitialMagnification','fit');
            set(gca,'Position',[0 0 1 1]);
            hold all;
        end    
        function plot_data(obj,x,y,varargin)
            [xi,yi] = obj.transform(x,y);
            plot(xi,yi,varargin{:});
        end
        function show_data_on_image(obj,x,y,varargin)
            % Create Figure
            obj.show_image;
                      
            % Plot Defined Points
            plot(obj.p1_cimg(1),obj.p1_cimg(2),'or','MarkerSize',5,'MarkerFaceColor','r')
            plot(obj.p2_cimg(1),obj.p2_cimg(2),'or','MarkerSize',5,'MarkerFaceColor','r')
            plot(obj.p3_cimg(1),obj.p3_cimg(2),'or','MarkerSize',5,'MarkerFaceColor','r')

            % Plot Origin
            obj.plot_data(0,0,'og','MarkerSize',5,'MarkerFaceColor','g')            
            
            % Plot Data
            if isempty(varargin)
                obj.plot_data(x,y,'b','LineWidth',2);
            else
                obj.plot_data(x,y,varargin{:});
            end
        end
    end
end

