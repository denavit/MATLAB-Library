classdef elasticAnalysis2d_member
    %elasticAnalysis2d_member Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        id                  % Unique integer identifies
        name = '';          % Optional friendly name in string format
        jointIDs            % ID numbers for I and J joints
        E                   % Modulus of Elasticity
        A                   % Cross Sectional Area
        I                   % Moment of Inertia
        endRelease = [0 0]; % Flag if moment is released or not (0 = not, 1 = released)
        numElements = 1;    % Number of elements used to model joint         
    end
    
    methods
        function obj=elasticAnalysis2d_member(id,jointIDs,E,A,I,endRelease)
            % elasticAnalysis2d_member(id,jointIDs,E,A,I,endRelease)
            obj.id = id;
            obj.jointIDs = jointIDs;
            obj.E = E;
            obj.A = A;
            obj.I = I;
            if nargin > 5
                obj.endRelease = endRelease;
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
        function obj = set.jointIDs(obj,jointIDs)
            if isequal(size(jointIDs),[2 1])
                obj.jointIDs = jointIDs';
            elseif isequal(size(jointIDs),[1 2])
                obj.jointIDs = jointIDs;
            else
                error('jointIDs is the wrong size')
            end            
        end
        function obj = set.E(obj,E)
            if isscalar(E) && isnumeric(E)
                obj.E = E;
            else
                error('E should be a numeric scalar');
            end
        end
        function obj = set.A(obj,A)
            if isscalar(A) && isnumeric(A)
                obj.A = A;
            else
                error('A should be a numeric scalar');
            end
        end
        function obj = set.I(obj,I)
            if isscalar(I) && isnumeric(I)
                obj.I = I;
            else
                error('I should be a numeric scalar');
            end
        end
        function obj = set.endRelease(obj,endRelease)
            % @todo - check that each element is either a zero or a one            
            if isequal(size(endRelease),[2 1])
                obj.endRelease = endRelease';
            elseif isequal(size(endRelease),[1 2])
                obj.endRelease = endRelease;
            else
                error('endRelease is the wrong size')
            end            
        end
        function obj = set.numElements(obj,numElements)
            if isscalar(numElements) && isnumeric(numElements) && round(numElements)==numElements
                obj.numElements = numElements;
            else
                error('numElements should be a scalar integer');
            end
        end        
        
        %% Output Function
        function print(obj,fid,flag)
            % print()
            % print(fid)
            % print(fid,flag)
            if (nargin < 2); fid = 1; end
            if (nargin < 3); flag = 1; end

            switch flag
                case 1
                    fprintf(fid,'    Member ID = %i, Joints (i,j) = (%i,%i), E = %g, A = %g, I = %g\n',...
                        obj.id,obj.jointIDs(1),obj.jointIDs(2),obj.E,obj.A,obj.I);
                otherwise 
                    warning('elasticAnalysis2d:badInput','Unknown print flag')
            end  
        end 
    end 
end

