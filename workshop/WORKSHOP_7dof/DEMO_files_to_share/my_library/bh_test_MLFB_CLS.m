classdef bh_test_MLFB_CLS
    
    properties
        TEST_MODEL = "bh_test_generic_MLFB_model"
    end
    
    methods
        function obj = bh_test_MLFB_CLS()
        end
    %======================================================================       
        function tf_is_config = isMexConfigured(obj)
                 tf_is_config = run_test_model(obj);
        end
    %======================================================================
    function tf_is_good = run_test_model(obj)
        
        tf_is_good  = false;
        
        fp_of_model = which("bh_test_generic_MLFB_model");

        [folder_of_model,~,~] = fileparts(fp_of_model);
        [folder_of_root,~,~]  = fileparts(folder_of_model);
         folder_of_test       =  [folder_of_root, filesep, 'Testing_sandpit'];
         
        my_pwd = pwd;
        
        if(isfolder(folder_of_test))
            cd(folder_of_test);
        else
            error("###_ERROR:  there is NO folder called <Testing_sandpit>");
        end
        
        if( isfolder("slprj") )
            rmdir("slprj","s");
        end
        
        % try and run the Simulink model that contains
        % a MATLAB Function block
        try
             y = sim( obj.TEST_MODEL );
             
             tf_is_good = true;
        catch
             tf_is_good = false;
        end   
        
        % back to where you were
        cd(my_pwd)
    end
    %======================================================================
    end % methods
end

