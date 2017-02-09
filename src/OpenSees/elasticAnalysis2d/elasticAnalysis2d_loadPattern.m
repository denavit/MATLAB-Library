classdef elasticAnalysis2d_loadPattern
    %elasticAnalysis2d_loadPattern Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        id              % Unique integer identifies
        name = '';      % Optional friendly name in string format
    end
    
    properties (SetAccess = private)
        jointLoads;
        memberLoads;
    end
    
    methods
        %% Constructor
        function obj=elasticAnalysis2d_loadPattern(id,name)
            % elasticAnalysis2d_loadPattern(id)
            % elasticAnalysis2d_loadPattern(id,name)
            obj.id = id;
            if ( nargin > 1 ); 
                obj.name = name;
            end
            obj.jointLoads.nodes = zeros(0,1);
            obj.jointLoads.loads = zeros(0,3);
            obj.memberLoads.members = zeros(0,1);
            obj.memberLoads.types = cell(0,1);
            obj.memberLoads.loads = cell(0,1);
        end
        
        %% Set Functions 
        function obj = set.id(obj,id)
            if isscalar(id)
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
        
        %% Information Functions
        function num = numJointLoads(obj)
            num = size(obj.jointLoads.nodes,1);
        end
        function num = numMemberLoads(obj)
            num = size(obj.memberLoads.members,1);
        end
        
        %% Adding Functions
        function newObj = addJointLoad(oldObj,nodes,loads)
            newObj = oldObj;
            if ( size(nodes,2) ~= 1 ); error('nodes should be a column vector'); end
            if ( size(loads,2) ~= 3 ); error('loads should have a width of 3'); end
            newObj.jointLoads.nodes = vertcat(newObj.jointLoads.nodes,nodes);
            newObj.jointLoads.loads = vertcat(newObj.jointLoads.loads,loads);
        end
        function newObj = addMemberLoad(oldObj,members,types,loads)
            newObj = oldObj;
            if ( size(members,2) ~= 1 ); error('members should be a column vector'); end
            if ischar(types); types = {types}; end;
            if ( size(types,2) ~= 1 ); error('types should be a column vector of cells'); end
            if isnumeric(loads)
                if size(loads,1)==1
                    loads = {loads};
                else
                    error('Invalid loads in addMemberLoad');
                end
            end
            if ( size(loads,2) ~= 1 ); error('loads should be a column vector of cells'); end
            newObj.memberLoads.members = vertcat(newObj.memberLoads.members,members);
            newObj.memberLoads.types = vertcat(newObj.memberLoads.types,types);
            newObj.memberLoads.loads = vertcat(newObj.memberLoads.loads,loads);
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
                    if strcmp(obj.name,'')
                        fprintf(fid,'    Load Pattern ID = %i\n',obj.id);
                    else
                        fprintf(fid,'    Load Pattern ID = %i, Name = %s\n',obj.id,obj.name);                        
                    end
                    if ~isempty(obj.jointLoads.nodes)
                        fprintf(fid,'      Number of Nodal Loads = %i\n',length(obj.jointLoads.nodes));
                        for i = 1:length(obj.jointLoads.nodes)
                            fprintf(fid,'        Node = %i, Load = (%g,%g,%g)\n',...
                                obj.jointLoads.nodes(i),obj.jointLoads.loads(i,1),...
                                obj.jointLoads.loads(i,2),obj.jointLoads.loads(i,3));
                        end
                    end
                    if ~isempty(obj.memberLoads.members)
                        fprintf(fid,'      Number of Element Loads = %i\n',length(obj.memberLoads.members));
                        for i = 1:length(obj.memberLoads.members)
                            fprintf(fid,'        Member = %i, Type = %s\n',...
                                obj.memberLoads.nodes(i),obj.memberLoads.type(i)); 
                        end
                    end                    
                otherwise 
                    warning('elasticAnalysis2d:badInput','Unknown print flag')
            end
        end       
    end
end

