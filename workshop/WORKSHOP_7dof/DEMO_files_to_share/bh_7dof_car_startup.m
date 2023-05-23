function bh_7dof_car_startup()
    
    % Add some folders to the MATLAB search path
    p                  = mfilename('fullpath');
    [folder,name,ext]  = fileparts(p);

    my_list_of_folders = { folder               , ...
                          'my_COMPONENT_systems', ...
                          'my_data'             , ...
                          'my_library'          , ...
                          'my_PICS'             , ...
                          'Task_road'           , ...
                          'Task_7dof_car'       , ...
                          };
    
    NUM_FOLDERS        = length(my_list_of_folders);
    
    for kk =1:NUM_FOLDERS
        if(1==kk)
            addpath(my_list_of_folders{kk},        '-begin');
        else
            addpath([folder, filesep, my_list_of_folders{kk}],'-begin');
        end      
    end

    % echo the first few elements of our search path
    sp             = path;
    TGT_SPLIT_CHAR = pathsep;
    C              = strsplit(sp, TGT_SPLIT_CHAR);
    
    fprintf('\n %s', repmat('*',1,50) );
    fprintf('\n Just added the following folders to the ');
    fprintf('\n HEAD of your search path: \n');
    fprintf('\n    ---> %s', C{1:NUM_FOLDERS});
    fprintf('\n %s', repmat('*',1,50) );
    fprintf('\n ... we are finished HERE ---> %s\n',mfilename);
        
    % assert that we have a new enough version to run this demo
    % R2020a is the minimum release and corresponds to MATLAB 
    % version 9.8
    MIN_required_ML_version = '9.8';
    if(verLessThan('MATLAB',MIN_required_ML_version))
          % inform the user that he needs a NEWER release and then exit
          error('###_ERROR:  you need at least R2019b to run this demo');
    end
       
    % OK, if you made it to here ... then Launch the 
    %     DEMO Navigation APP
    bh_7dof_car_content_navigator_app
    
end % function
