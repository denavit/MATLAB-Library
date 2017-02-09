function modelFileData = elasticAnalysis2d_writeModelInputFile(obj,fid,options)

%% Retreive what is needed from input
% Required input
loadFactors         = options.loadFactors;
analysisType        = options.analysis_type;
scratch_path        = options.scratch_path;

% Optional input
if isfield(options,'notional_load_ratio')
    notionalLoadRatio = options.notional_load_ratio;
else
    notionalLoadRatio = 0;
end
if isfield(options,'stiffness_reduction')
    stiffnessReduction = options.stiffness_reduction;
else
    stiffnessReduction = 1.0;
end    

%% Initilize node and member counters
nNode = 1;
nElement = 1;
nodeForJoint = zeros(obj.numJoints,1);
forMember = struct;

%% Open file and write header and basic model info
fprintf(fid, '# elasticAnalysis2d, Model: %s\n\n',obj.name);
fprintf(fid, 'model BasicBuilder -ndm 2 -ndf 3 \n\n');

%% Nodes and Boundary Conditions
fprintf(fid, '# Nodes and Boundary Conditions \n');
for i = 1:obj.numJoints
    fprintf(fid, 'node %i %g %g \n',nNode,obj.joints(i).coords(1),obj.joints(i).coords(2));
    if ~isequal(obj.joints(i).bc,[0 0 0])
        fprintf(fid, 'fix %i %i %i %i \n',nNode,obj.joints(i).bc(1),obj.joints(i).bc(2),obj.joints(i).bc(3));
    end
    nodeForJoint(i) = nNode; % Record the node number used for the joint
    nNode = nNode+1;
end
fprintf(fid, '\n');

%% Equal DOF Constraints
fprintf(fid, '# Equal DOF Constraints \n');
for i = 1:obj.numEqualDofs
    rNode = nodeForJoint(obj.findJoint(obj.equalDofs(i).retainedJoint));
    cNode = nodeForJoint(obj.findJoint(obj.equalDofs(i).constrainedJoint));
    dofs = find(obj.equalDofs(i).dofs==1);
    fprintf(fid, 'equalDOF %i %i %s \n',rNode,cNode,num2str(dofs));
end
fprintf(fid, '\n');

%% Geometric Transformation
fprintf(fid, '# Geometric Transformation \n');
fprintf(fid, 'set transfTag 1\n');
switch lower(analysisType)
    case {'linear','firstorder'}
        fprintf(fid, 'geomTransf Linear $transfTag \n\n');
    case {'nonlinear','secondorder'}
        fprintf(fid, 'geomTransf PDelta $transfTag \n\n');
    case 'exact'
        fprintf(fid, 'geomTransf Corotational $transfTag \n\n');        
    otherwise
        error('Unknown analysisType');
end

%% Elements
fprintf(fid, '# Elements \n');
for i = 1:obj.numMembers
    iJoint = obj.findJoint(obj.members(i).jointIDs(1));
    jJoint = obj.findJoint(obj.members(i).jointIDs(2));
    iCoords = obj.joints(iJoint).coords;
    jCoords = obj.joints(jJoint).coords;
    numElements = obj.members(i).numElements;
    iE = obj.members(i).E*stiffnessReduction;
    iA = obj.members(i).A;
    iI = obj.members(i).I;
    if (numElements == 1)
        % Create node for i end if joint i is released and add equal dof contraint
        if ( obj.members(i).endRelease(1) == 1 )
             iNode = nNode;
             nNode = nNode+1;
             fprintf(fid, 'node %i %g %g \n',iNode,iCoords);
             fprintf(fid, 'equalDOF %i %i 1 2 \n',nodeForJoint(iJoint),iNode);
        else
             iNode = nodeForJoint(iJoint);
        end
        % Create node for j end if joint j is released and add equal dof contraint
        if ( obj.members(i).endRelease(2) == 1 )
             jNode = nNode;
             nNode = nNode+1;
             fprintf(fid, 'node %i %g %g \n',jNode,jCoords);             
             fprintf(fid, 'equalDOF %i %i 1 2 \n',nodeForJoint(jJoint),jNode);
        else
             jNode = nodeForJoint(jJoint);
        end        
        fprintf(fid, 'element elasticBeamColumn %i %i %i %g %g %g $transfTag \n',...
            nElement,iNode,jNode,iA,iE,iI);
        forMember(i).node = [iNode jNode];
        forMember(i).element = nElement;
        nElement = nElement+1;          
    else
        % Find Nodal Coords
        nodeCoords = horzcat(linspace(iCoords(1),jCoords(1),numElements+1)',...
            linspace(iCoords(2),jCoords(2),numElements+1)');
        % First Element
        % Create node for i end if joint i is released and add equal dof contraint
        if ( obj.members(i).endRelease(1) == 1 )
             iNode = nNode;
             nNode = nNode+1;
             fprintf(fid, 'node %i %g %g \n',iNode,iCoords);
             fprintf(fid, 'equalDOF %i %i 1 2 \n',nodeForJoint(iJoint),iNode);
        else
             iNode = nodeForJoint(iJoint);
        end   
        % Create Node for j end
        jNode = nNode;
        nNode = nNode+1;
        fprintf(fid, 'node %i %g %g \n',jNode,nodeCoords(2,:));
        % Create element
        fprintf(fid, 'element elasticBeamColumn %i %i %i %g %g %g $transfTag \n',...
            nElement,iNode,jNode,iA,iE,iI);
        forMember(i).node = [iNode jNode];
        forMember(i).element = nElement;
        nElement = nElement+1;
        
        % Middle Elements
        for j = 2:numElements-1
            iNode = jNode;
            % Create Node for j end
            jNode = nNode;
            nNode = nNode+1;
            fprintf(fid, 'node %i %g %g \n',jNode,nodeCoords(j+1,:));
            % Create element
            fprintf(fid, 'element elasticBeamColumn %i %i %i %g %g %g $transfTag \n',...
                nElement,iNode,jNode,iA,iE,iI);
            forMember(i).node = horzcat(forMember(i).node,jNode);
            forMember(i).element = horzcat(forMember(i).element,nElement);
            nElement = nElement+1;
        end
        
        % Last Element
        iNode = jNode;
        % Create node for j end if joint j is released and add equal dof contraint
        if ( obj.members(i).endRelease(2) == 1 )
             jNode = nNode;
             nNode = nNode+1;
             fprintf(fid, 'node %i %g %g \n',jNode,jCoords);             
             fprintf(fid, 'equalDOF %i %i 1 2 \n',nodeForJoint(jJoint),jNode);
        else
             jNode = nodeForJoint(jJoint);
        end  
        fprintf(fid, 'element elasticBeamColumn %i %i %i %g %g %g $transfTag \n',...
            nElement,iNode,jNode,iA,iE,iI);
        forMember(i).node = horzcat(forMember(i).node,jNode);
        forMember(i).element = horzcat(forMember(i).element,nElement);
        nElement = nElement+1;
    end
end
fprintf(fid, '\n');

%% Loads
fprintf(fid, '# Loads \n');
for i = 1:size(loadFactors,1)
    ilpID = loadFactors(i,1);
    ilp   = obj.findLoadPattern(ilpID);
    fprintf(fid, 'timeSeries Linear %i -factor %g\n',ilp,loadFactors(i,2));
    fprintf(fid, 'pattern Plain %i %i { \n',ilp,ilp);
    % Nodal Loads
    for j = 1:length(obj.loadPatterns(ilp).jointLoads.nodes)
        Fx = obj.loadPatterns(ilp).jointLoads.loads(j,1);
        Fy = obj.loadPatterns(ilp).jointLoads.loads(j,2);
        Fq = obj.loadPatterns(ilp).jointLoads.loads(j,3);
        fprintf(fid, '  load %i %g %g %g \n',...
            nodeForJoint(obj.findJoint(obj.loadPatterns(ilp).jointLoads.nodes(j))),...
            Fx-notionalLoadRatio*Fy,Fy,Fq);
    end
    % Element Load
    for j = 1:size(obj.loadPatterns(ilp).memberLoads.types,1)
        member = obj.loadPatterns(ilp).memberLoads.members(j);
        memberID = obj.findMember(member);
        numElements = obj.members(memberID).numElements; 
        type   = obj.loadPatterns(ilp).memberLoads.types{j};
        loads  = obj.loadPatterns(ilp).memberLoads.loads{j};

        [iJointCoords,jJointCoords] = obj.memberCoords(member);
        angle = atan2(jJointCoords(2)-iJointCoords(2),jJointCoords(1)-iJointCoords(1));
        
        switch type
            case 'Uniform'
                eles = num2str(forMember(memberID).element);
                if numel(loads) == 1
                    Wzo = loads;
                    Wxo = 0.0;
                elseif numel(loads) == 2
                    Wzo = loads(1);
                    Wxo = loads(2);
                else
                    error('Invalid size for loads in Uniform Load');
                end
                notionalLoad = -notionalLoadRatio*(cos(angle)*Wzo+sin(angle)*Wxo);
                Wz = Wzo - sin(angle)*notionalLoad;
                Wx = Wxo + cos(angle)*notionalLoad;
                fprintf(fid, '  eleLoad -ele %s -type -beamUniform %g %g \n',...
                    eles,Wz,Wx);
            case 'Point'
                if numel(loads) == 2
                    Pzo = loads(1);
                    Pxo = 0.0;
                    xL = loads(2);
                elseif numel(loads) == 3
                    Pzo = loads(1);
                    Pxo = loads(2);
                    xL = loads(3);
                else
                    error('Invalid size for loads in Point Load');
                end
                notionalLoad = -notionalLoadRatio*(cos(angle)*Pzo+sin(angle)*Pxo);
                Pz = Pzo - sin(angle)*notionalLoad;
                Px = Pxo + cos(angle)*notionalLoad;
                xL = xL*numElements;
                xLi = floor(xL);
                xL = xL-xLi;
                if (xL == 0)
                    iJoint = obj.findJoint(obj.members(memberID).jointIDs(1));
                    jJoint = obj.findJoint(obj.members(memberID).jointIDs(2));
                    iJointCoords = obj.joints(iJoint).coords;
                    jJointCoords = obj.joints(jJoint).coords;
                    VEC = jJointCoords-iJointCoords;
                    [Q,~] = cart2pol(VEC(1),VEC(2));
                    T = [cos(Q) -sin(Q); sin(Q) cos(Q)];
                    P = T*[Px Pz]';
                    fprintf(fid, '  load %i %g %g 0.0 \n',...
                        forMember(memberID).node(xLi+1),P(1),P(2));
                else
                    fprintf(fid, '  eleLoad -ele %i -type -beamPoint %g %g %g \n',...
                        forMember(memberID).element(xLi+1),Pz,xL,Px);                    
                end
            otherwise
                error('Unknown member load type');
        end
    end
    fprintf(fid, '}\n');
end
fprintf(fid, '\n');

%% Recorders
fprintf(fid, '# Recorders \n');
% Joint Displacements and Reactions (One file for all elements)
recorderFilenames.jointDisplacements = fullfile(scratch_path,'elasticAnalysis2d_jointDisplacements.out');
recorderFilenames.jointReactions     = fullfile(scratch_path,'elasticAnalysis2d_jointReactions.out');
nodes = sprintf('%i ',nodeForJoint);
fprintf(fid, 'recorder Node -file %s -node %s -dof 1 2 3 disp\n',...
    path_for_tcl(recorderFilenames.jointDisplacements),nodes);
fprintf(fid, 'recorder Node -file %s -node %s -dof 1 2 3 reaction\n',...
    path_for_tcl(recorderFilenames.jointReactions),nodes);

% Member Displacement and Forces (One file for each member)
for i = 1:obj.numMembers
    % Local forces for each member
    recorderFilenames.memberForces{i} = fullfile(scratch_path,sprintf('elasticAnalysis2d_memberForces%i.out',i));
    elements = sprintf('%i ',forMember(i).element);
    fprintf(fid, 'recorder Element -file {%s} -ele %s localForces \n',...
        path_for_tcl(recorderFilenames.memberForces{i}),elements);
    
    % Displacements at the nodes for each member
    recorderFilenames.memberDisplacements{i} = fullfile(scratch_path,sprintf('elasticAnalysis2d_memberDisplacements%i.out',i));
    nodes = sprintf('%i ',forMember(i).node);
    fprintf(fid, 'recorder Node    -file {%s} -node %s -dof 1 2 3 disp \n',...
        path_for_tcl(recorderFilenames.memberDisplacements{i}),nodes);
end
fprintf(fid, '\n');

% Store data from input file. 
modelFileData.nodeForJoint = nodeForJoint;
modelFileData.forMember = forMember;
modelFileData.recorderFilenames = recorderFilenames;

end



function tclpath = path_for_tcl(path)
tclpath = regexprep(path,'\','/');
end