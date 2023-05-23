classdef bh_7dof_car_demo_selector_CLS
   
    properties
        ROOT_FOLDER = '';
        REF_FOLDER  = '';
    end
    
methods
function OBJ = bh_7dof_car_demo_selector_CLS()
    p = mfilename('fullpath');
    [pathstr,~,~] = fileparts(p);
    % am assuming THIS file is IN a folder under the ROOT folder !
    [pathstr,~,~]   = fileparts(pathstr);
    OBJ.ROOT_FOLDER = string(pathstr);
    OBJ.REF_FOLDER  = string(pathstr);
end
%------------------------------------------------------------------
function clear_base(OBJ)
    LOC_clear_base();
end
%------------------------------------------------------------------
function close_editor(OBJ)
    LOC_close_M_editor();
end
%------------------------------------------------------------------
function open_content_app(OBJ)   
    %bh_3LINK_IK_content_navigator_app();   
end
%_#########################################################################
%   S  Y  S  T  E  M       C  H  E  C  K
%_#########################################################################

function do_system_check(obj, app)
    LOC_do_system_check(app);
end
%------------------------------------------------------------------
%_#########################################################################
%   R O A D 
%_#########################################################################
function open_ROAD_demos(OBJ, DEMO_ID)
    DEMO_ID      = upper(string(DEMO_ID));
    
    % all should be UPPER case
    allowed_string_list = [ "NOTHING", "ROAD_HUMPS_EXPLORE"];
        
    % confirm our selection is a member of our allowed list
    mustBeMember(DEMO_ID, allowed_string_list)

    % make sure you're in the WORKING folder
    OBJ.cd_into_working_folder();
    
    % cd into the Local task folder
    cd("Task_road")
    
    % clear BASE workspace
    LOC_clear_base();  
   
    switch DEMO_ID
        case "NOTHING"
              return
        case "ROAD_HUMPS_EXPLORE"
             edit('DEMO_bh_test_roads_profile.mlx');
    end % switch
end
%------------------------------------------------------------------
%_#################################################################
%   7 - D O F      C A R
%_#################################################################
function open_7DOFCAR_demos(OBJ, DEMO_ID)
    DEMO_ID      = upper(string(DEMO_ID));
    
    % all should be UPPER case
    allowed_string_list = [ "NOTHING", ...
                            "DERIVE_EOMS", ...
                            "BRAIN_MODEL", ...
                            "SIMSCAPE_MODEL", ...
                            "SIMULATE_AT_MODAL", ...
                            "SIMULATE_ROAD_HUMPS", ...
                            "SIMULATE_VALIDATE"];
        
    % confirm our selection is a member of our allowed list
    mustBeMember(DEMO_ID, allowed_string_list)

    % make sure you're in the WORKING folder
    OBJ.cd_into_working_folder();
    
    % cd into the Local task folder
    cd("Task_7dof_car")
    
    % clear BASE workspace
    LOC_clear_base();  
   
    switch DEMO_ID
        case "NOTHING"
              return
              
        case "DERIVE_EOMS"
             edit('DEMO_01_bh_7dof_derivation.mlx');
             
        case "BRAIN_MODEL"     
             open_system('bh_COMPONENT_model_7dof_car_HAND.slx');  
            
        case "SIMSCAPE_MODEL"    
             open_system('bh_COMPONENT_model_7dof_car_SCAPE.slx'); 
            
        case "SIMULATE_AT_MODAL"
             edit('DEMO_02_bh_7dof_launch_sim_modal.mlx');
             
        case "SIMULATE_ROAD_HUMPS"
             edit('DEMO_03_bh_7dof_launch_sim_humps.mlx');
             
        case "SIMULATE_VALIDATE"
            edit("bh_VALIDATE_7dof_car_HAND_vs_SIMSCAPE.mlx");
            
    end % switch
end
%_#################################################################
%   D A T A
%_#################################################################
function open_DATA_demos(OBJ, DEMO_ID)
    DEMO_ID      = upper(string(DEMO_ID));
    
    % all should be UPPER case
    allowed_string_list = ["EXCEL_DATABASE_FILE"];
        
    % confirm our selection is a member of our allowed list
    mustBeMember(DEMO_ID, allowed_string_list)

    % make sure you're in the WORKING folder
    OBJ.cd_into_working_folder();
       
    switch DEMO_ID
        case "EXCEL_DATABASE_FILE"
              winopen("bh_7dof_car_database.xlsx")
    end % switch
end
%_#################################################################
%   A C T I V E    S U S P E N S I O N
%_#################################################################
function open_active_suspension_demos(OBJ, DEMO_ID)
    DEMO_ID      = upper(string(DEMO_ID));
    
    % all should be UPPER case
    allowed_string_list = [ "NOTHING", ...
                            "DESIGN"  ];
        
    % confirm our selection is a member of our allowed list
    mustBeMember(DEMO_ID, allowed_string_list)

    % make sure you're in the WORKING folder
    OBJ.cd_into_working_folder();
    
    % cd into the Local task folder
    cd("Task_ActiveSuspension")
    
    % clear BASE workspace
    LOC_clear_base();  
   
    switch DEMO_ID
        case "NOTHING"
              return
              
        case "DESIGN"
             edit('DEMO_Active_Suspension.mlx');           
    end % switch
end


%_#################################################################
%   T I P S  and  F Y I
%_#################################################################
function open_tips_fyi_topic(OBJ, DEMO_ID)
    DEMO_ID      = upper(string(DEMO_ID));
    
    % all should be UPPER case
    allowed_string_list = ["SUBSYSTEM_REFERENCE_01"];
        
    % confirm our selection is a member of our allowed list
    mustBeMember(DEMO_ID, allowed_string_list)
    
    switch DEMO_ID
        case "SUBSYSTEM_REFERENCE_01"
              edit("bh_tip_01_modelling_SUBSYSTEM_REFERENCE.mlx"); 
    end % switch
    
end
%_#################################################################
%   U T I L I T I E S
%_#################################################################
function wk_folder = get_working_folder(OBJ)   
    wk_folder = string(OBJ.ROOT_FOLDER);
end
%------------------------------------------------------------------
function cd_into_working_folder(OBJ)
    my_pwf                  = string(pwd);   
    THE_TGT_WORKING_FOLDER  = OBJ.ROOT_FOLDER;
    
    if(my_pwf ~= THE_TGT_WORKING_FOLDER)
         cd( THE_TGT_WORKING_FOLDER );
    end
end
%------------------------------------------------------------------

end % methods
%_*************************************************************************
end % classdef
%_#########################################################################
function LOC_clear_base()
   evalin('base', 'clear all; clc');
   disp(' ...### BASE workspace has been cleared')
end
%**************************************************************************
function LOC_close_M_editor
%close the M-file Editor
  evalin('base',['com.mathworks.mlservices.MatlabDesktopServices.getDesktop().closeGroup(''Editor'')']);
end
%**************************************************************************
function LOC_do_system_check(app)
    DEFAULT_COLOR        = [255 191 0]/255; % orange
    app.LAMP_ML.Color    = DEFAULT_COLOR;
    app.LAMP_SYMT.Color  = DEFAULT_COLOR;            
    app.LAMP_SL.Color    = DEFAULT_COLOR;            
    app.LAMP_SC.Color    = DEFAULT_COLOR;            
    app.LAMP_SCMB.Color  = DEFAULT_COLOR;

    % just pause for a second
    pause(1);

    % what products do you have
    v_S          = ver;
    prod_list_CE = {v_S.Name}';
    %prod_list_CE([31, 63, 67]) = []; % for testing only !
    % update the lamps
    my_database_CE = ...
        {
          'MATLAB',                    'LAMP_ML';
          'Symbolic Math Toolbox',     'LAMP_SYMT';   
          'Simulink',                  'LAMP_SL';            
          'Simscape',                  'LAMP_SC';
          'Simscape Multibody',        'LAMP_SCMB';
          };
      NUM_LAMPS = size(my_database_CE,1);
      for kk=1:NUM_LAMPS
          the_prod_str  = my_database_CE{kk,1};
          the_field_str = my_database_CE{kk,2};
          tf_i_have_prod = strcmp(prod_list_CE, the_prod_str);

          if(1==nnz(tf_i_have_prod))
              app.(the_field_str).Color = 'green';
          else
              app.(the_field_str).Color = 'red';
          end
      end % for
      
      % What about having an installed C-compiler ?
      % - we need this for the models containing MATLAb Function BLOCKS
            
      app.LAMP_C_Comp_A.Color = DEFAULT_COLOR;
      app.LAMP_C_Comp_B.Color = DEFAULT_COLOR;
      app.LAMP_C_Comp_C.Color = DEFAULT_COLOR;
      
      pause(2)
      
      app.LAMP_C_Comp_A.Color = 'green'; pause(2)
      app.LAMP_C_Comp_B.Color = 'green'; pause(2)
           
      test_OBJ = bh_test_MLFB_CLS;
      
      if(test_OBJ.isMexConfigured)
          app.LAMP_C_Comp_C.Color = 'green';
      else
          app.LAMP_C_Comp_C.Color = 'red';
      end
      
end % function LOC_do_system_check(app)    end
%**************************************************************************

