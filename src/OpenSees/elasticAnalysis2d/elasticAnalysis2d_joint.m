classdef elasticAnalysis2d_joint
    %elasticAnalysis2d_joint Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        id              % Unique integer identifies
        name = '';      % Optional friendly name in string format
        coords          % Joint coordinates in the two dimensions
        bc = [0 0 0];   % Flag if dof is fixed or free (0 = free, 1 = fixed)
    end
    
    methods
        %% Constructor
        function obj=elasticAnalysis2d_joint(id,coords,bc)
            % elasticAnalysis2d_joint(id,coords,bc)
            obj.id = id;
            obj.coords = coords; 
            if nargin > 2
                obj.bc = bc;
            end
        end
        
        %% Set Functions
        function obj = set.id(obj,id)
            if isscalar(id) && isnumeric(id)
                obj.id = id;
            else
                error('id should be a numeric scalar');
            end
        end
        function obj = set.name(obj,name)
            if ischar(name)
                obj.name = name;
            else
                error('name should be a character string');
            end
        end
        function obj = set.coords(obj,coords)
            if isequal(size(coords),[2 1])
                obj.coords = coords';
            elseif isequal(size(coords),[1 2])
                obj.coords = coords;
            else
                error('coords is the wrong size')
            end            
        end
        function obj = set.bc(obj,bc)
            % @todo - check that each element is either a zero or a one
            if isequal(size(bc),[3 1])
                obj.bc = bc';
            elseif isequal(size(bc),[1 3])
                obj.bc = bc;
            else
                error('bc is the wrong size')
            end            
        end 
        
        %% Output Functions
        function print(obj,fid,flag)
            % print()
            % print(fid)
            % print(fid,flag)
            if (nargin < 2); fid = 1; end
            if (nargin < 3); flag = 1; end

            switch flag
                case 1
                    fprintf(fid,'    Joint ID: %i, Coords: (%g,%g)\n',...
                        obj.id,obj.coords(1),obj.coords(2));
                otherwise 
                    warning('elasticAnalysis2d:badInput','Unknown print flag')
            end
        end 
    end
    
end

