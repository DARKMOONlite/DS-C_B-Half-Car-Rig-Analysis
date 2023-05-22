
% directory = 'LAB TEST FOLDER\Centre\Damped\'
% x = lvm_import("LAB TEST FOLDER\Centre\Damped\test 2_2, damped 500.lvm");
% 
% img_dir = dir(fullfile(directory,'*.lvm'))
% 

classdef DataReader

    properties 
        filelist
    end
    methods 
        function obj = DataReader(directory)
            obj.filelist = dir(fullfile(directory,'*.lvm'));
    
        end

    end
    
    
end