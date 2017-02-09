classdef elasticAnalysis2d_model < handle
    %ELASTICANALYSIS2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name                      % Friendly name in string format
        defaultNumElements = [];   % Default numElements
    end
    
    properties (SetAccess=private)
        joints                        % Vector of joint objects
        jointIDs = zeros(0,1);        % Vector containg IDs of the joints
        equalDofs                     % Vector of equalDof objects       
        members                       % Vector of member objects
        memberIDs = zeros(0,1);       % Vector containg IDs of the members
        loadPatterns                  % Vector of loadPattern objects
        loadPatternIDs = zeros(0,1);  % Vector containg IDs of the load patterns
    end    
    
    properties (Hidden = true)
        backgroundColor = [0.75 0.75 0.75];
        passColor = [0.0 1.0 0.0];
        failColor = [1.0 0.0 0.0];
    end
    
    methods
        %% Constructer
        function obj=elasticAnalysis2d_model(name,defaultNumElements)
            % elasticAnalysis2d_joint(name)
            obj.name = name;
            if nargin > 1
                obj.defaultNumElements = defaultNumElements;
            end
        end
        
        
        %% Set Functions
        function set.name(obj,name)
            if ischar(name)
                obj.name = name;
            else
                error('name should be a character string');
            end
        end   
        function set.defaultNumElements(obj,defaultNumElements)
            if isscalar(defaultNumElements)  && isnumeric(defaultNumElements) ...
                    && round(defaultNumElements)==defaultNumElements
                obj.defaultNumElements = defaultNumElements;
            else
                error('defaultNumElements should be a scalar integer');
            end
        end
        
        %% Add and Remove Functions
        function addJoint(obj,varargin)
            newJoint = elasticAnalysis2d_joint(varargin{:});
            newId = newJoint.id;
            if ~isempty(obj.findJoint(newId))
                error('New Joint ID is not unique');
            end
            obj.joints = vertcat(obj.joints,newJoint);
            obj.jointIDs = vertcat(obj.jointIDs,newId);
        end
        function removeJoint(obj,jointID)
            jointNum = obj.findJoint(jointID);
            assert(~isempty(jointNum),'Joint ID %g not defined',jointID)
            obj.joints = vertcat(obj.joints(1:(jointNum-1)),obj.joints((jointNum+1):end));
            obj.jointIDs = vertcat(obj.jointIDs(1:(jointNum-1)),obj.jointIDs((jointNum+1):end));
        end
        function addEqualDof(obj,varargin)
            obj.equalDofs = vertcat(obj.equalDofs,elasticAnalysis2d_equalDof(varargin{:}));
        end        
        function addMember(obj,varargin)
            newMember = elasticAnalysis2d_member(varargin{:});
            if ~isempty(obj.defaultNumElements)
                newMember.numElements = obj.defaultNumElements;
            end
            newId = newMember.id;
            if ~isempty(obj.findMember(newId))
                error('New Member ID is not unique');
            end            
            obj.members = vertcat(obj.members,newMember);
            obj.memberIDs = vertcat(obj.memberIDs,newId);            
        end
        function removeMember(obj,memberID)
            memberNum = obj.findMember(memberID);
            assert(~isempty(memberNum),'Member ID %g not defined',memberID)
            obj.members = vertcat(obj.members(1:(memberNum-1)),obj.members((memberNum+1):end));
            obj.memberIDs = vertcat(obj.memberIDs(1:(memberNum-1)),obj.memberIDs((memberNum+1):end));
        end        
        function addLoadPattern(obj,varargin)
            newLoadPattern = elasticAnalysis2d_loadPattern(varargin{:});
            newId = newLoadPattern.id;
            if ~isempty(obj.findLoadPattern(newId))
                error('New LoadPattern ID is not unique');
            end            
            obj.loadPatterns = vertcat(obj.loadPatterns,newLoadPattern);
            obj.loadPatternIDs = vertcat(obj.loadPatternIDs,newId); 
        end    
        
        %% Modify Properties
        function setJointProperty(obj,jointID,propName,varargin)
            jointNum = obj.findJoint(jointID);
            switch propName
                case 'id'
                    obj.joints(jointNum).id = varargin{1};
                case 'name'
                    obj.joints(jointNum).name = varargin{1};
                case 'coords'
                    obj.joints(jointNum).coords = varargin{1};
                case 'bc'
                    obj.joints(jointNum).bc = varargin{1};
                otherwise
                    error('Unknown property name');
            end
        end
        function x = getJointProperty(obj,jointID,propName)
            jointNum = obj.findJoint(jointID);
            switch propName
                case 'id'
                    x = obj.joints(jointNum).id;
                case 'name'
                    x = obj.joints(jointNum).name;
                case 'coords'
                    x = obj.joints(jointNum).coords;
                case 'bc'
                    x = obj.joints(jointNum).bc;
                otherwise
                    error('Unknown property name');
            end
        end        
        function setMemberProperty(obj,member,propName,varargin)
            memberID = obj.findMember(member);
            switch propName
                case 'id'
                    obj.members(memberID).id = varargin{1};
                case 'name'
                    obj.members(memberID).name = varargin{1};
                case 'jointIDs'
                    obj.members(memberID).jointIDs = varargin{1};
                case 'section'
                    obj.members(memberID).section = varargin{1};
                case 'endRelease'
                    obj.members(memberID).endRelease = varargin{1};
                case 'outOfStraightness'
                    obj.members(memberID) = ...
                        obj.members(memberID).setOutOfStraightness(varargin{:});  
                case 'numElements'
                    obj.members(memberID).numElements = varargin{1};                    
                otherwise
                    error('Unknown property name');
            end
        end
        function addJointLoad(obj,loadPattern,nodes,loads)
            iLP = obj.findLoadPattern(loadPattern);
            obj.loadPatterns(iLP) = obj.loadPatterns(iLP).addJointLoad(nodes,loads);
        end     
        function addMemberLoad(obj,loadPattern,members,types,loads)
            iLP = obj.findLoadPattern(loadPattern);
            obj.loadPatterns(iLP) = obj.loadPatterns(iLP).addMemberLoad(members,types,loads);
        end        
        
        %% Information Functions
        function num = numJoints(obj)
            num = length(obj.joints);
        end
        function num = numEqualDofs(obj)
            num = length(obj.equalDofs);
        end   
        function num = numMembers(obj)
            num = length(obj.members);
        end        
        function num = numLoadPatterns(obj)
            num = length(obj.loadPatterns);
        end
        function num = findJoint(obj,id)
            if numel(id) == 1
                num = find(obj.jointIDs==id); 
            else
                num = zeros(size(id));
                for i = 1:numel(id)
                    n = find(obj.jointIDs==id(i));
                    if isempty(n)
                        error('Invalid joint id');
                    else
                        num(i) = n;
                    end
                end
            end
        end
        function num = findMember(obj,id) 
            num = find(obj.memberIDs==id);          
        end    
        function num = findLoadPattern(obj,id)
            num = find(obj.loadPatternIDs==id);           
        end   
        function coords = jointCoords(obj,id)
            num = obj.findJoint(id);
            coords = obj.joints(num).coords;
        end
        function [iJointCoords,jJointCoords] = memberCoords(obj,memberID)
            memberNum = obj.findMember(memberID);
            iJoint = obj.findJoint(obj.members(memberNum).jointIDs(1));
            jJoint = obj.findJoint(obj.members(memberNum).jointIDs(2));
            iJointCoords = obj.joints(iJoint).coords;
            jJointCoords = obj.joints(jJoint).coords;
        end
        function [xMin,xMax,yMin,yMax] = structureBounds(obj)
            xMin =  Inf;
            xMax = -Inf;
            yMin =  Inf;
            yMax = -Inf;
            for i = 1:obj.numJoints
                coords = obj.joints(i).coords;
                if ( coords(1)<xMin ); xMin = coords(1); end
                if ( coords(1)>xMax ); xMax = coords(1); end
                if ( coords(2)<yMin ); yMin = coords(2); end
                if ( coords(2)>yMax ); yMax = coords(2); end
            end
        end
        
        
        %% Output Functions
        function print(obj,fid,flag)
            % print()
            % print(fid)
            % print(fid,flag)
            if (nargin < 2); fid = 1; end
            if (nargin < 3); flag = 1; end

            fprintf(fid,'elasticAnalysis2d Model:\n');
            fprintf(fid,'  Name = %s\n',obj.name);
            fprintf(fid,'  Number of Joints = %i\n',obj.numJoints);
            for i = 1:obj.numJoints
                obj.joints(i).print(fid,flag);
            end
            fprintf(fid,'  Number of Equal DOF Contraints = %i\n',obj.numEqualDofs);     
            fprintf(fid,'  Number of Members = %i\n',obj.numMembers);
            for i = 1:obj.numMembers
                obj.members(i).print(fid,flag);
            end           
            fprintf(fid,'  Number of Load Patterns = %i\n',obj.numLoadPatterns);
            for i = 1:obj.numLoadPatterns
                obj.loadPatterns(i).print(fid,flag);
            end               
        end
        
        %% Plotting Functions
        function plotJoints(obj,color)
            if nargin < 2
                color = obj.backgroundColor;
            end
            if strcmpi(color,'backrgound')
                color = obj.backgroundColor;
            end
            for i = 1:length(obj.joints)
                plot(obj.joints(i).coords(1),obj.joints(i).coords(2),'o',...
                    'Color',color,'MarkerFaceColor',color,...
                    'MarkerSize',7)
            end
        end
        function plotMembers(obj,varargin)
            structuralModel2d_plotMembers(obj,varargin{:});
        end
        function plotEndReleases(obj,color)
            if nargin < 2
                color = obj.backgroundColor;
            end
            if strcmpi(color,'backrgound')
                color = obj.backgroundColor;
            end
            for i = 1:obj.numMembers
                iJoint = obj.findJoint(obj.members(i).jointIDs(1));
                jJoint = obj.findJoint(obj.members(i).jointIDs(2));
                iJointCoords = obj.joints(iJoint).coords;
                jJointCoords = obj.joints(jJoint).coords;
                if ( obj.members(i).endRelease(1) == 1 )
                    coords = iJointCoords + 0.1*(jJointCoords-iJointCoords);
                    plot(coords(1),coords(2),'o','LineWidth',2,...
                        'Color',color)
                end
                if ( obj.members(i).endRelease(2) == 1 )
                    coords = iJointCoords + 0.9*(jJointCoords-iJointCoords);
                    plot(coords(1),coords(2),'o','LineWidth',2,...
                        'Color',color)
                end
            end
        end        
        function plotMemberEnds(obj)
            for i = 1:obj.numMembers
                iJoint = obj.findJoint(obj.members(i).jointIDs(1));
                jJoint = obj.findJoint(obj.members(i).jointIDs(2));
                iJointCoords = obj.joints(iJoint).coords;
                jJointCoords = obj.joints(jJoint).coords;
                coords = iJointCoords + 0.1*(jJointCoords-iJointCoords);
                h1 = plot(coords(1),coords(2),'or','LineWidth',2);
                coords = iJointCoords + 0.9*(jJointCoords-iJointCoords);
                h2 = plot(coords(1),coords(2),'ob','linewidth',2);
            end
            legend([h1 h2],{'I End','J End'})
        end        
        function plotLoadPattern(obj,loadPatternID)
            structuralModel2d_plotLoadPattern(obj,loadPatternID);
        end
        function plotJointIDs(obj)
            for i = 1:length(obj.joints)
                text(obj.joints(i).coords(1),obj.joints(i).coords(2),...
                    sprintf('%i',obj.joints(i).id),...
                    'VerticalAlignment','middle',...
                    'HorizontalAlignment','center');
            end           
        end
        function plotConstraints(obj,color)
            if nargin < 2
                color = obj.backgroundColor;
            end
            if strcmpi(color,'backrgound')
                color = obj.backgroundColor;
            end
            [xMin,xMax,yMin,yMax] = obj.structureBounds;
            w = 0.1*max([xMax-xMin yMax-yMin]);
            for i = 1:length(obj.joints)
                ibc = obj.joints(i).bc;
                icoords = obj.joints(i).coords;
                if isequal(ibc,[0 0 0])
                    continue;
                elseif isequal(ibc,[1 0 0])
                    [x,y] = shapes.rollerSupport(icoords(1),icoords(2),0,w);
                elseif isequal(ibc,[0 1 0])
                    [x,y] = shapes.rollerSupport(icoords(1),icoords(2),-pi/2,w);
                elseif isequal(ibc,[1 1 0])
                    [x,y] = shapes.pinSupport(icoords(1),icoords(2),-pi/2,w);
                elseif isequal(ibc,[1 1 1])
                    [x,y] = shapes.rigidSupport(icoords(1),icoords(2),-pi/2,w);
                else
                    warning('Constraint plot not implemented');
                end
                plot(x,y,'LineWidth',2,'Color',color)
                    
            end            
        end        
        function plotMemberIDs(obj)
            obj.plotMemberData(obj.memberIDs);         
        end     
        function plotMemberData(obj,memberData,format)
            if nargin < 3
                format = [];
            end
            for i = 1:obj.numMembers
                iJoint = obj.findJoint(obj.members(i).jointIDs(1));
                jJoint = obj.findJoint(obj.members(i).jointIDs(2));
                iJointCoords = obj.joints(iJoint).coords;
                jJointCoords = obj.joints(jJoint).coords;
                coords = (jJointCoords+iJointCoords)/2;
                
                if iscell(memberData)
                    if isempty(memberData{i})
                        iText = '';
                    elseif ~isempty(format)
                        iText = sprintf(format,memberData{i});
                    elseif isinteger(memberData{i})
                        iText = sprintf('%i',memberData{i});
                    elseif isnumeric(memberData{i})
                        iText = sprintf('%g',memberData{i});
                    elseif ischar(memberData{i})
                        iText = memberData{i};
                    else
                        error('Unknown type for memberData');
                    end
                elseif isnumeric(memberData)
                    if isnan(memberData(i))
                        iText = '';
                    elseif ~isempty(format)
                        iText = sprintf(format,memberData(i));
                    else
                        iText = num2str(memberData(i));
                    end
                else
                    error('Unknown type for memberData');
                end                   
                text(coords(1),coords(2),iText,...
                    'VerticalAlignment','middle',...
                    'HorizontalAlignment','center');                
            end
        end        
        function plotStructure(obj)
            figure
            hold all
            obj.plotJoints('backrgound');
            if nargin < 2
                obj.plotMembers;
            else
                if ischar(memberValues)
                    memberValues = obj.getMemberData(memberValues);
                end
                if isnumeric(memberValues)
                    obj.plotMembers('ColorValues',memberValues);
                    colorbar
                else
                    obj.plotMembers;
                end
                if nargin < 3
                    obj.plotMemberData(memberValues);
                else
                    obj.plotMemberData(memberValues,format);
                end
                
            end            
            obj.setAxisLimits(0.1);
        end
        function setAxisLimits(obj,margin)
            [xMin,xMax,yMin,yMax] = obj.structureBounds;
            xlim([xMin-margin*(xMax-xMin) xMax+margin*(xMax-xMin)])
            ylim([yMin-margin*(yMax-yMin) yMax+margin*(yMax-yMin)])            
        end
        
        function modelFileData = writeModelInputFile(obj,fid,options)
            modelFileData = elasticAnalysis2d_writeModelInputFile(obj,fid,options);
        end
    end
    
end

