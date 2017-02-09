classdef OpenSeesAnalysis < handle
    
    properties
        pathOfOpenSees
        echoOpenSeesOutput = false;
        deleteFilesAfterAnalysis = true;
        scratchPath;
    end
    
    methods
        % Constructor
        function obj = OpenSeesAnalysis
            try
                obj.pathOfOpenSees = pathOf.OpenSees;
            catch
                obj.pathOfOpenSees = '';
            end
            try
                obj.scratchPath = pathOf.scratch;
            catch
                obj.scratchPath = '';
            end
        end
        
        % Set and Get Functions
        function set.pathOfOpenSees(obj,pathOfOpenSees)
            if ischar(pathOfOpenSees)
                obj.pathOfOpenSees = pathOfOpenSees;
            else
                error('pathOfOpenSees should be a character string');
            end
        end
        function pathOfOpenSees = get.pathOfOpenSees(obj)
            if isempty(obj.pathOfOpenSees)
                error('pathOfOpenSees not set');
            end
            pathOfOpenSees = obj.pathOfOpenSees;
        end
        function set.scratchPath(obj,scratchPath)
            if ischar(scratchPath)
                obj.scratchPath = scratchPath;
            else
                error('scratchPath should be a character string');
            end
        end
        
        % Analysis Functions
        function str = scratchFile(obj,filename)
            str = fullfile(obj.scratchPath,filename);
        end
        function [status, result] = runOpenSees(obj,inputFile,echo)
            if nargin < 3
                echo = obj.echoOpenSeesOutput;
            end
            command = sprintf('"%s" "%s"',obj.pathOfOpenSees,inputFile);
            if echo
                [status, result] = system(command,'-echo');
            else
                [status, result] = system(command);
            end
        end
        
    end
    
    methods (Static)
        function tclpath = path_for_tcl(path)
            tclpath = regexprep(path,'\','/');
        end
    end
end