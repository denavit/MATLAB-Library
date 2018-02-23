classdef SectionAnalysis < OpenSeesAnalysis
    
    properties
        section_def
        section_tag = 1;
    end
    
    methods
        %% Constructor
        function obj = SectionAnalysis(def,tag)
            obj.section_def = def;
            if nargin > 1
                obj.section_tag = tag;
            end
        end
        
        %% Set and get functions
        function set.section_def(obj,section_def)
            if ischar(section_def)
                if exist(section_def,'file') == 2
                    fid = fopen(section_def);
                    A = textscan(fid,'%s\n','Whitespace','\n');
                    obj.section_def = A{1};
                else
                    obj.section_def = {section_def};
                end
            elseif iscell(section_def) && isvector(section_def)
                obj.section_def = section_def;
            else
                error('section_def should be a cell vector');
            end
        end
        
        %% Section Discretization Functions
        function [fiberLocZ,fiberLocY,fiberArea,fiberMat] = getDiscretization(obj)
            
            % Filenames to use
            filename_input = obj.scratchFile('SectionAnalyis_getDiscretization_InputFile.tcl');
            filename_print = obj.scratchFile('SectionAnalyis_getDiscretization_PrintFile.txt');
            
            % Create Tcl File
            myFile = fopen(filename_input,'w');
            fprintf(myFile, 'model BasicBuilder -ndm 3 -ndf 6 \n');
            for i = 1:length(obj.section_def)
                fprintf(myFile, '%s\n',obj.section_def{i});
            end
            fprintf(myFile, 'node 1 0.0 0.0 0.0\n');
            fprintf(myFile, 'node 2 0.0 0.0 0.0\n');            
            fprintf(myFile, 'element zeroLengthSection 1 1 2 1 \n');
            fprintf(myFile, 'print {%s} -ele -flag 3 1\n',obj.path_for_tcl(filename_print));
            fprintf(myFile, 'exit 2 \n');
            fclose(myFile);
                       
            % Run OpenSees
            if exist(filename_print,'file') == 2
                delete(filename_print);
            end
            [status, result] = obj.runOpenSees(filename_input);
            switch status
                case 2
                    % analysis successful
                otherwise
                    fprintf('%s\n',result);
                    error('Analysis Failed in Unknown Manner (exit code: %i)',status);
            end            

            % Read Results           
            C = dlmread(filename_print,' ',4,0);
            fiberMat  = C(:,1);
            fiberLocY = C(:,2);
            fiberLocZ = C(:,3);
            fiberArea = C(:,4);
            
            % Clean Folder
            if obj.deleteFilesAfterAnalysis
                delete(filename_input,filename_print);
            end
        end     
        function plotDiscretization(obj,varargin)
            % Set default values
            plotAs2dSection = false;
            
            % Read arguments
            for i = 1:(nargin-1)
                switch varargin{i}
                    case {'PlotAs2d','plotAs2dSection'}
                        plotAs2dSection  = true;  
                    case {'PlotAs3d','plotAs3dSection'}
                        plotAs2dSection  = false;                  
                    otherwise
                end
                
            end            

            % Run Analysis
            [fiberLocZ,fiberLocY,fiberArea,fiberMat] = obj.getDiscretization; 
            numMats = length(unique(fiberMat));
            if numMats <= 7
                colormap(lines(numMats))
            else
                temp = hsv(numMats);
                colormap(vertcat(temp(1:2:end,:),temp(2:2:end,:)));
            end
            if plotAs2dSection
                fiberLocZ = 1:length(fiberLocY);
            end
            scatter(fiberLocZ,fiberLocY,20,fiberMat,'filled')
        end
        function printDiscretization(obj,fileID)
            if nargin < 2
                fileID = 1;
            end
            
            % Run Analysis
            [fiberLocZ,fiberLocY,fiberArea,fiberMat] = obj.getDiscretization;
            
            % Print Data
            fprintf(fileID,'z location\ty location\tarea\tmaterial id\n');
            for i = 1:length(fiberArea)
                fprintf(fileID,'%g\t%g\t%g\t%i\n',...
                    fiberLocZ(i),fiberLocY(i),fiberArea(i),fiberMat(i));
            end
        end
        function printMaterialInfo(obj,fileID)
            if nargin < 2
                fileID = 1;
            end
            
            % Run Analysis
            [fiberLocZ,fiberLocY,fiberArea,fiberMat] = obj.getDiscretization;        
            
            % Print data
            fprintf(fileID,'  Material  |  # Fibers  |    Area    |     Iz     |     Iy     \n');
            fprintf(fileID,'------------+------------+------------+------------+------------\n');        
            uniqueFiberMat = unique(fiberMat);
            for i = 1:length(uniqueFiberMat)
                ind = find( fiberMat == uniqueFiberMat(i) );
                PartSectionArea = sum(fiberArea(ind));
                PartSectionIy   = sum(fiberArea(ind).*fiberLocZ(ind).^2);
                PartSectionIz   = sum(fiberArea(ind).*fiberLocY(ind).^2);
                fprintf(fileID,'     %-7.i|   %-9.i|% -12.6g|% -12.6g|% -12.6g\n',...
                    uniqueFiberMat(i),numel(ind),PartSectionArea,PartSectionIz,PartSectionIy);
            end            
            fprintf(fileID,'------------+------------+------------+------------+------------\n'); 
            SectionArea = sum(fiberArea);
            SectionIy   = sum(fiberArea.*fiberLocZ.^2);
            SectionIz   = sum(fiberArea.*fiberLocY.^2);
            fprintf(fileID,'    Total   |   %-9.i|% -12.6g|% -12.6g|% -12.6g\n',...
                numel(fiberArea),SectionArea,SectionIz,SectionIy);
        end
        
        %% 2d Interaction Functions
        function results = runSectionToPeak2d(obj,loadingType,deformationStep,maxNumSteps,P,numStepsAxial)

            % Filenames
            filename_input        = fullfile(obj.scratchPath,'sectionAnalysis_Input.tcl');
            filename_output_force = fullfile(obj.scratchPath,'sectionAnalysis_OutputForce.out');
            filename_output_def   = fullfile(obj.scratchPath,'sectionAnalysis_OutputDisp.out');
            filename_output_eigen = fullfile(obj.scratchPath,'sectionAnalysis_OutputEigen.out');

            % Create .tcl file
            fid = fopen(filename_input,'w');
            fprintf(fid,'model BasicBuilder -ndm 2 -ndf 3 \n');
            fprintf(fid,'package require OpenSeesComposite \n');
            fprintf(fid,'namespace import OpenSeesComposite::* \n');

            % Node and boundry conditions
            fprintf(fid,'node 1 0.0 0.0 \n');
            fprintf(fid,'node 2 0.0 0.0 \n');
            fprintf(fid,'mass 1 1.0 1.0 1.0 \n');
            fprintf(fid,'mass 2 1.0 1.0 1.0 \n');
            fprintf(fid,'fix 1 1 1 1 \n');
            fprintf(fid,'fix 2 0 1 0 \n');

            % Define the Section
            if iscell(obj.section_def)
                for i = 1:length(obj.section_def)
                    fprintf(fid,'%s \n',obj.section_def{i});
                end
            else
                fprintf(fid,'%s \n',obj.section_def);
            end
            fprintf(fid,'element zeroLengthSection 1 1 2 1 \n');

            % Recorders
            fprintf(fid,'recorder Node -file {%s} -node 1 -dof 1 2 3 reaction \n',filename_output_force);
            fprintf(fid,'recorder Node -file {%s} -node 2 -dof 1 2 3 disp \n',filename_output_def);

            % One Step with Nothing
            fprintf(fid,'system UmfPack \n');
            fprintf(fid,'constraints Transformation \n');
            fprintf(fid,'test NormDispIncr 1.0e-6 10 1 \n');
            fprintf(fid,'algorithm Newton \n');
            fprintf(fid,'numberer RCM \n');
            fprintf(fid,'integrator LoadControl 0\n');
            fprintf(fid,'analysis Static\n');
            fprintf(fid,'set ok [analyze 1]\n');
            fprintf(fid,'if {$ok == 0} {\n');
            fprintf(fid,'  set lowestEigen [eigenRecorder {%s} 1 -generalized -fullGenLapack] \n',filename_output_eigen);
            fprintf(fid,'}\n');
            
            switch loadingType
                case 'AxialOnly'
                    % Axial Loading
                    fprintf(fid,'timeSeries Linear 1 \n');
                    fprintf(fid,'pattern Plain 1 1 { \n');
                    fprintf(fid,'  load 2 1.0 0.0 0.0 \n');
                    fprintf(fid,'} \n');
                    fprintf(fid,'integrator DisplacementControl 2 1 %g \n',deformationStep);
                    fprintf(fid,'analysis Static \n');
                    fprintf(fid,'for {set i 1} {$i <= %i} {incr i} { \n',maxNumSteps);
                    fprintf(fid,'  set ok [analyze 1] \n');
                    fprintf(fid,'  if {$ok != 0} { exit 3 } \n');
                    fprintf(fid,'  set lowestEigen [eigenRecorder {%s} 1 -generalized -fullGenLapack] \n',filename_output_eigen);
                    fprintf(fid,'  if {$lowestEigen < 0.0} { exit 1 } \n');
                    fprintf(fid,'} \n');
                    fprintf(fid,'exit 2 \n');
                case 'NonProportional'
                    % Axial Loading
                    fprintf(fid,'timeSeries Linear 1 \n');
                    fprintf(fid,'pattern Plain 1 1 { \n');
                    fprintf(fid,'  load 2 1.0 0.0 0.0 \n');
                    fprintf(fid,'} \n');
                    fprintf(fid,'integrator LoadControl %g \n',P/numStepsAxial);
                    fprintf(fid,'analysis Static \n');
                    fprintf(fid,'for {set i 1} {$i <= %i} {incr i} { \n',numStepsAxial);
                    fprintf(fid,'  set ok [analyze 1] \n');
                    fprintf(fid,'  if {$ok != 0} { exit 4 } \n');
                    fprintf(fid,'  set lowestEigen [eigenRecorder {%s} 1 -generalized -fullGenLapack] \n',filename_output_eigen);
                    fprintf(fid,'  if {$lowestEigen < 0.0} { exit 5 } \n');
                    fprintf(fid,'} \n');
                    fprintf(fid,'loadConst -time 0.0 \n');
                    
                    % Bending Loading
                    fprintf(fid,'timeSeries Linear 2 \n');
                    fprintf(fid,'pattern Plain 2 2 { \n');
                    fprintf(fid,'  load 2 0.0 0.0 1.0 \n');
                    fprintf(fid,'} \n');
                    fprintf(fid,'integrator DisplacementControl 2 3 %g \n',deformationStep);
                    fprintf(fid,'analysis Static \n');
                    fprintf(fid,'for {set i 1} {$i <= %i} {incr i} { \n',maxNumSteps);
                    fprintf(fid,'  set ok [analyze 1] \n');
                    fprintf(fid,'  if {$ok != 0} { exit 3 } \n');
                    fprintf(fid,'  set lowestEigen [eigenRecorder {%s} 1 -generalized -fullGenLapack] \n',filename_output_eigen);
                    fprintf(fid,'  if {$lowestEigen < 0.0} { exit 1 } \n');
                    fprintf(fid,'} \n');  
                    fprintf(fid,'exit 2 \n');
                otherwise
                    error('Unknown loadingType: %s',loadingType);
            end

            fclose(fid);

            % Run OpenSees
            if exist(filename_output_eigen,'file')
                delete(filename_output_eigen)
            end
            [status, result] = obj.runOpenSees(filename_input);
            
            results = struct;
            results.textOutput = result;
            
            switch status
                case 1
                    % Analysis successful
                    results.status = 'PeakObtained';
                case 2
                    % Deformation limits reached
                    results.status = 'DeformationLimits';
                case 3
                    results.status = 'AnalysisFailed';
                    warning('sectionInspector:AnalysisResult','Analysis failed in main loading');
                case 4
                    fprintf('%s\n',result);
                    error('Analysis failed in axial loading');
                case 5
                    fprintf('%s\n',result);
                    error('Peak reached during axial loading');
                otherwise
                    fprintf('%s\n',result);
                    error('Analysis ended in an unknown manner, exit code = %i',status);
            end

            % Read Results
            disp  = dlmread(filename_output_def);
            force = -1*dlmread(filename_output_force);
            eigen  = dlmread(filename_output_eigen);
            switch loadingType
                case 'AxialOnly'
                    startPoint = 1;
                case 'NonProportional'
                    startPoint = numStepsAxial;
                otherwise
                    error('Unknown loadingType');
            end     
            results.axialStrain = disp(startPoint:end,1);
            results.curvature   = disp(startPoint:end,3);
            results.axialForce  = force(startPoint:end,1);
            results.moment      = force(startPoint:end,3);
            results.eigen       = eigen(startPoint:end,2);
            
            % Clean Folder
            if obj.deleteFilesAfterAnalysis
                delete(filename_input,filename_output_force,filename_output_def,...
                    filename_output_eigen);
            end
        end
        function [P,M,status] = sectionInteraction2d(obj,numPoints,...
                deformationStep_Axial,maxNumSteps_Axial,...
                deformationStep_Bending,maxNumSteps_Bending)
            
            % Initilize results
            M = zeros(numPoints,1);
            P = zeros(numPoints,1);
            status = cell(numPoints,1);
            
            % Run axial only analysis to get Pn
            results = obj.runSectionToPeak2d(deformationStep_Axial,...
                maxNumSteps_Axial);
            M(1) = results.limitPoint_M;
            P(1) = results.limitPoint_P;
            status{1} = results.status;
            
            % Run nonproportional analyses at different axial loads
            Ps = linspace(P(1),0,numPoints);
            for i = 2:length(Ps)
                results = obj.runSectionToPeak2d(...
                    deformationStep_Bending,maxNumSteps_Bending,Ps(i),10);
                M(i) = results.limitPoint_M;
                P(i) = results.limitPoint_P;
                status{i} = results.status;
            end
            
            if nargout < 3
                clear status;
            end
        end
        
        %% One Dimensional Stress Strain Analysis
        function results = getStressStrain(obj,peakPoints,matID,rateType,rateValue)
            %
            % getStressStrain(peakPoints)
            % getStressStrain(peakPoints,matID)
            % getStressStrain(peakPoints,matID,rateType,rateValue)
            %
            %
            
            if nargin < 3
                matID = 1;
            end
            
            if nargin < 4
                rateType  = 'Steps';
                rateValue = 100;
            end
            
            switch rateType
                case 'StrainRate'
                    rate = rateValue;
                case 'Steps'
                    rate = sum(abs(diff(peakPoints)))/rateValue;
                otherwise
                    error('Unknown rateType');
            end          
            
            % Filenames To Use
            % (these files will be overwriten and then deleted)
            filename_input =        fullfile(obj.scratchPath,'runMaterialAnalysisInputFile.tcl');
            filename_disp =         fullfile(obj.scratchPath,'runMaterialAnalysisPattern.txt');
            filename_output_force = fullfile(obj.scratchPath,'runMaterialAnalysisPatternForce.out');
            filename_output_def =   fullfile(obj.scratchPath,'runMaterialAnalysisPatternDisp.out');
            filename_output_stiff = fullfile(obj.scratchPath,'runMaterialAnalysisPatternStiff.out');
            
            % Generate imposed displacment history
            numbers = fillOutNumbers(peakPoints,rate);
            numbers = vertcat(numbers(1),numbers,numbers(end));
            numSteps = length(numbers)-2;
            save(filename_disp,'numbers','-ASCII');
            
            % Create .tcl file
            myFile = fopen(filename_input,'w');
            fprintf(myFile, 'model BasicBuilder -ndm 2 -ndf 3 \n');
            fprintf(myFile, 'node 1 0.0 0.0 \n');
            fprintf(myFile, 'node 2 1.0 0.0 \n');
            fprintf(myFile, 'fix 1 1 1 1 \n');
            fprintf(myFile, 'fix 2 0 1 1 \n');
            for i = 1:length(obj.section_def)
                fprintf(myFile, '%s\n',obj.section_def{i});
            end
            fprintf(myFile, 'element truss 1 1 2 1.0 %i \n',matID);
            fprintf(myFile, 'timeSeries Path 1 -dt 1.0 -filePath %s \n',OpenSeesAnalysis.path_for_tcl(filename_disp));
            fprintf(myFile, 'pattern Plain 1 1 { \n');
            fprintf(myFile, '    sp 2 1 1.0 \n');
            fprintf(myFile, '} \n');
            fprintf(myFile, 'recorder Element -file %s -ele 1  force \n',OpenSeesAnalysis.path_for_tcl(filename_output_force));
            fprintf(myFile, 'recorder Node    -file %s -node 2 -dof 1 disp \n',OpenSeesAnalysis.path_for_tcl(filename_output_def));
            fprintf(myFile, 'recorder Element -file %s -ele 1  stiff \n',OpenSeesAnalysis.path_for_tcl(filename_output_stiff));
            fprintf(myFile, 'system UmfPack \n');
            fprintf(myFile, 'constraints Penalty 1.0e12 1.0e12 \n');
            fprintf(myFile, 'test NormDispIncr 1.0e-8 10 0 \n');
            fprintf(myFile, 'algorithm Newton \n');
            fprintf(myFile, 'numberer RCM \n');
            fprintf(myFile, 'integrator LoadControl 1.0 \n');
            fprintf(myFile, 'analysis Static \n');
            fprintf(myFile, ['analyze ' num2str(numSteps) ' \n']);
            fprintf(myFile, 'exit 1 \n');
            fclose(myFile);
            
            % Run OpenSees
            [status, result] = obj.runOpenSees(filename_input);
            switch status
                case 1
                    % Analysis Successful
                otherwise
                    fprintf('%s\n',result);
            end
            
            
            % Read Results
            results.disp = dlmread(filename_output_def);
            temp = dlmread(filename_output_force);
            results.force = temp(:,4);
            try
                results.stiff = dlmread(filename_output_stiff);
            catch
                results.stiff = [];
            end
            
            % Clean Folder
            delete(filename_input,filename_disp,filename_output_force,...
                filename_output_def,filename_output_stiff);
        end
        
        %% Section Analysis
        function results = runSectionAnalysis(obj,dofs,peakPoints,rate,varargin)
            %
            % dofs is a vector defining whether each dof is:
            %     0 - not tested (displacement = 0)
            %     1 - displacement control
            %     2 - force control
            %     3 - natural boundary conditions (force = 0)
            % the length of the vector is:
            %     3 (axial,shear,bending) for two-dimensional analyses
            %     6 (axial,shearY,shearZ,torsion,bendingY,bendingZ) for
            %          three-dimensional analyses
            %
            
            % Set default values
            getFiberData = false;
            
            % Read arguments
            for i = 1:(nargin-4)
                switch lower(varargin{i})
                    case 'getfiberdata'
                        getFiberData = true;  
                    otherwise
                        error('Unknown argument: %s',varargin{i});
                end
            end             
            
            % Check input and set infomation variables
            if isequal(size(dofs),[3 1]) || isequal(size(dofs),[1 3])
                ndm = 2;
                ndf = 3;
            elseif isequal(size(dofs),[6 1]) || isequal(size(dofs),[1 6])
                ndm = 3;
                ndf = 6;
            else
                error('dofs is the wrong size');
            end
                        
            dofFixed    = find(dofs==0);
            dofDisp     = find(dofs==1);
            dofForce    = find(dofs==2);
            dofNatural  = find(dofs==3);
            numdofControlled = length(dofDisp)+length(dofForce);
            if (length(dofFixed)+length(dofDisp)+length(dofForce)+length(dofNatural) ~= length(dofs))
                error('dofs is not formatted correctly');
            end
            
            if ( size(peakPoints,2) ~= numdofControlled )
                error('peakPoints is the wrong size');
            end
            
            if ( size(rate,2) ~= numdofControlled )
                error('rate is the wrong size');
            end
            
            % Filenames
            filename_input        = obj.scratchFile('runSectionAnalysis_Input.tcl');
            filename_output_force = obj.scratchFile('runSectionAnalysis_OutputForce.out');
            filename_output_def   = obj.scratchFile('runSectionAnalysis_OutputDisp.out');
            filepattern_pattern   = 'runSectionAnalysis_Pattern_%i.txt';
            filepattern_fiber     = 'runSectionAnalysis_Fiber_%i.txt';
            
            % Retreive and Store Fiber Data
            if getFiberData
                [fiberLocZ,fiberLocY,fiberArea,fiberMat] = obj.getDiscretization;
                numFibers = length(fiberLocZ);
                % Store Fiber Data
                results.fiber_data.num_fibers = numFibers;
                results.fiber_data.z_pos  = fiberLocZ;
                results.fiber_data.y_pos  = fiberLocY;
                results.fiber_data.area   = fiberArea;
                results.fiber_data.mat_id = fiberMat;
            end            
            
            % Generate imposed displacment history
            numbers = fillOutNumbers(peakPoints,rate);
            numbers = vertcat(numbers(1,:),numbers,numbers(end,:));
            numSteps = size(numbers,1)-2;
            j = 1;
            for i = find( dofs==1 | dofs==2 )
                filename = obj.scratchFile(sprintf(filepattern_pattern,i));
                iNumbers = numbers(:,j);
                j = j + 1;
                save(filename,'iNumbers','-ASCII');
            end
            
            % Create .tcl file
            myFile = fopen(filename_input,'w');
            
            % Node and boundry conditions
            if (ndm == 2)
                fprintf(myFile, 'model BasicBuilder -ndm 2 -ndf 3 \n');
                fprintf(myFile, 'node 1 0.0 0.0 \n');
                fprintf(myFile, 'node 2 0.0 0.0 \n');
                fprintf(myFile, 'fix 1 1 1 1 \n');
                fprintf(myFile, 'fix 2 %i %i %i \n',dofs==0);
            elseif (ndm == 3)
                fprintf(myFile, 'model BasicBuilder -ndm 3 -ndf 6 \n');
                fprintf(myFile, 'node 1 0.0 0.0 0.0 \n');
                fprintf(myFile, 'node 2 0.0 0.0 0.0 \n');
                fprintf(myFile, 'fix 1 1 1 1 1 1 1 \n');
                fprintf(myFile, 'fix 2 %i %i %i %i %i %i \n',dofs==0);
            end
            
            % Define the Section
            for i = 1:length(obj.section_def)
                fprintf(myFile, '%s\n',obj.section_def{i});
            end
            fprintf(myFile, 'element zeroLengthSection 1 1 2 1 \n');
            
            % Loading
            for i = 1:ndf
                switch dofs(i)
                    case 1
                        ifilename = obj.scratchFile(sprintf(filepattern_pattern,i));
                        %fprintf(myFile,'set pattern_path "%s" \n',obj.path_for_tcl(ifilename));
                        fprintf(myFile,'pattern Plain %i "Series -dt 1.0 -filePath {%s} -factor 1.0" { \n',i,obj.path_for_tcl(ifilename));
                        fprintf(myFile,'    sp 2 %i 1.0 \n',i);
                        fprintf(myFile,'} \n');
                    case 2
                        ifilename = obj.scratchFile(sprintf(filepattern_pattern,i));
                        fprintf(myFile,'pattern Plain %i "Series -dt 1.0 -filePath {%s} -factor 1.0" { \n',i,obj.path_for_tcl(ifilename));
                        if ( ndm == 2 )
                            temp = zeros(1,3);
                            temp(i) = 1;
                            fprintf(myFile, '    load 2 %i %i %i \n',temp);
                        elseif ( ndm == 3 )
                            temp = zeros(1,6);
                            temp(i) = 1;
                            fprintf(myFile, '    load 2 %i %i %i %i %i %i \n',temp);
                        end
                        fprintf(myFile, '} \n');
                    otherwise
                        
                end
            end
            
            % Recorders
            if (ndm == 2)
                fprintf(myFile,'recorder Node -file {%s} -node 1 -dof 1 2 3 reaction \n',obj.path_for_tcl(filename_output_force));
                fprintf(myFile,'recorder Node -file {%s} -node 2 -dof 1 2 3 disp \n',obj.path_for_tcl(filename_output_def));
            elseif (ndm == 3)
                fprintf(myFile,'recorder Node -file {%s} -node 1 -dof 1 2 3 4 5 6 reaction \n',obj.path_for_tcl(filename_output_force));
                fprintf(myFile,'recorder Node -file {%s} -node 2 -dof 1 2 3 4 5 6 disp \n',obj.path_for_tcl(filename_output_def));
            end
            if getFiberData
                for i = 1:numFibers
                    ifilename = obj.scratchFile(sprintf(filepattern_fiber,i));
                    fprintf(myFile,'recorder Element -file "%s" -ele 1 section fiber %i stressANDstrain \n',obj.path_for_tcl(ifilename),i-1);
                    %fprintf(myFile,'recorder Element -file "%s" -ele 1 section fiber %g %g %i stressANDstrain \n',obj.path_for_tcl(ifilename),...
                    %    results.fiber_data.y_pos(i),results.fiber_data.z_pos(i),results.fiber_data.mat_id(i));
                end
            end
            
            % Analysis
            fprintf(myFile, 'system UmfPack \n');
            fprintf(myFile, 'constraints Penalty 1.0e12 1.0e12 \n');
            fprintf(myFile, 'test NormDispIncr 1.0e-8 10 0 \n');
            fprintf(myFile, 'algorithm Newton \n');
            fprintf(myFile, 'numberer RCM \n');
            fprintf(myFile, 'integrator LoadControl 1.0 \n');
            fprintf(myFile, 'analysis Static \n');
            fprintf(myFile, 'set ok [analyze %i] \n',numSteps);
            fprintf(myFile, 'if {$ok != 0} {exit 2} \n');
            fprintf(myFile, 'exit 1 \n');
            fclose(myFile);
            
            % Run OpenSees
            [status, result] = obj.runOpenSees(filename_input);
            results.textOutput  = result;
                        
            switch status
                case 1
                    results.status = 'Analysis Successful';
                case 2
                    results.status = 'Analysis Failed';
                    fprintf('%s\n',result);
                    warning('SectionInspector:Analysis','Analysis failed\n');
                otherwise
                    fprintf('%s\n',result);
                    error('Analysis ended in an unknown manner, exit code = %i',status);
            end
                        
            % Read Results
            results.disp  = dlmread(filename_output_def,' ');
            results.force = -1*dlmread(filename_output_force,' ');
            if (ndm == 2)
                results.axialStrain = results.disp(:,1);
                results.shearStrain = results.disp(:,2);
                results.curvature   = results.disp(:,3);
                results.axialForce  = results.force(:,1);
                results.shearForce  = results.force(:,2);
                results.moment      = results.force(:,3);
            elseif (ndm == 3)
                results.axialStrain  = results.disp(:,1);
                results.shearYStrain = results.disp(:,2);
                results.shearZStrain = results.disp(:,3);
                results.twist        = results.disp(:,4);
                results.curvatureY   = results.disp(:,5);
                results.curvatureZ   = results.disp(:,6);
                results.axialForce   = results.force(:,1);
                results.shearYForce  = results.force(:,2);
                results.shearZForce  = results.force(:,3);
                results.torsion      = results.force(:,4);
                results.momentY      = results.force(:,5);
                results.momentZ      = results.force(:,6);
            end
            if getFiberData
                results.fiber_results(length(results.fiber_data.z_pos)) = struct;
                for i = 1:numFibers
                    ifilename = obj.scratchFile(sprintf(filepattern_fiber,i));
                    C = dlmread(ifilename,' ');
                    results.fiber_results(i).stress = C(:,1);
                    results.fiber_results(i).strain = C(:,2);
                end
            end
            
            % Clean Folder
            if obj.deleteFilesAfterAnalysis
                delete(filename_input,filename_output_force,filename_output_def);
                for i = find( dofs==1 | dofs==2 )
                    ifilename = obj.scratchFile(sprintf(filepattern_pattern,i));
                    delete(ifilename)
                end
                if getFiberData
                    for i = 1:numFibers
                        ifilename = obj.scratchFile(sprintf(filepattern_fiber,i));
                        delete(ifilename)
                    end
                end
            end
        end
    end
end
