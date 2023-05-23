classdef bh_7dof_car_database_CLS
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = protected)
        Filename
        Sheetname_params
        Sheetname_IC
        car_tab
        IC_tab
        %------------------------------------------------------------------
        m_c
        I_cZ
        I_cX
        m_bR
        m_fR
        m_bL
        m_fL
        %------------------------------------------------------------------
        k_bR
        k_fR
        k_bL
        k_fL
        c_bR
        c_fR
        c_bL
        c_fL
        %------------------------------------------------------------------
        k_bRG
        k_fRG
        k_bLG
        k_fLG
        c_bRG
        c_fRG
        c_bLG
        c_fLG
        %------------------------------------------------------------------
        L_wheelbase
        L_track
        %------------------------------------------------------------------
        L_f
        L_b
        L_R
        L_L
        %------------------------------------------------------------------
        % initial conditions
        %------------------------------------------------------------------
        IC_pos
        IC_vel

    end
%==========================================================================    
    methods
        function OBJ = bh_7dof_car_database_CLS(filename, sheetname_pars, sheetname_ICs)
            arguments
                filename       (1,1) string
                sheetname_pars (1,1) string
                sheetname_ICs  (1,1) string
            end
            
            OBJ.Filename  = filename;
            
            % read the car PARAMETERS  
            %  --> configure table such that Name can be used as an index mechanism
            OBJ.Sheetname_params = sheetname_pars;
            OBJ.car_tab          = readtable(filename, "Sheet", sheetname_pars, "TextType", "string");           
            OBJ.car_tab.Properties.RowNames = OBJ.car_tab.Name;  
                
            OBJ = LOC_assign_car_params(OBJ);
            
            % read the INITIAL CONDITIONS
            %  --> configure table such that Name can be used as an index mechanism
            OBJ.Sheetname_IC = sheetname_ICs;
            OBJ.IC_tab       = readtable(filename, "Sheet", sheetname_ICs, "TextType", "string");
            OBJ.IC_tab.Properties.RowNames = OBJ.IC_tab.Name;  
            
            OBJ = LOC_assign_ICs(OBJ);
        end
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        %------------------------------------------------------------------
    end
end
%_#########################################################################
%_ LOCAL SUBFUNCTIONS
%_#########################################################################
function OBJ = LOC_assign_car_params(OBJ, tab)
    OBJ.m_c         = OBJ.car_tab{"m_c"        , "Value"};
    OBJ.I_cZ        = OBJ.car_tab{"I_cZ"       , "Value"};
    OBJ.I_cX        = OBJ.car_tab{"I_cX"       , "Value"};
    OBJ.m_bR        = OBJ.car_tab{"m_bR"       , "Value"};
    OBJ.m_fR        = OBJ.car_tab{"m_fR"       , "Value"};
    OBJ.m_bL        = OBJ.car_tab{"m_bL"       , "Value"};
    OBJ.m_fL        = OBJ.car_tab{"m_fL"       , "Value"};
    OBJ.k_bR        = OBJ.car_tab{"k_bR"       , "Value"};
    OBJ.k_fR        = OBJ.car_tab{"k_fR"       , "Value"};
    OBJ.k_bL        = OBJ.car_tab{"k_bL"       , "Value"};
    OBJ.k_fL        = OBJ.car_tab{"k_fL"       , "Value"};
    OBJ.c_bR        = OBJ.car_tab{"c_bR"       , "Value"};
    OBJ.c_fR        = OBJ.car_tab{"c_fR"       , "Value"};
    OBJ.c_bL        = OBJ.car_tab{"c_bL"       , "Value"};
    OBJ.c_fL        = OBJ.car_tab{"c_fL"       , "Value"};
    OBJ.k_bRG       = OBJ.car_tab{"k_bRG"      , "Value"};
    OBJ.k_fRG       = OBJ.car_tab{"k_fRG"      , "Value"};
    OBJ.k_bLG       = OBJ.car_tab{"k_bLG"      , "Value"};
    OBJ.k_fLG       = OBJ.car_tab{"k_fLG"      , "Value"};
    OBJ.c_bRG       = OBJ.car_tab{"c_bRG"      , "Value"};
    OBJ.c_fRG       = OBJ.car_tab{"c_fRG"      , "Value"};
    OBJ.c_bLG       = OBJ.car_tab{"c_bLG"      , "Value"};
    OBJ.c_fLG       = OBJ.car_tab{"c_fLG"      , "Value"};
    OBJ.L_wheelbase = OBJ.car_tab{"L_wheelbase", "Value"};
    OBJ.L_track     = OBJ.car_tab{"L_track"    , "Value"};
    OBJ.L_f         = OBJ.car_tab{"L_f"        , "Value"};
    OBJ.L_b         = OBJ.car_tab{"L_b"        , "Value"};
    OBJ.L_R         = OBJ.car_tab{"L_R"        , "Value"};
    OBJ.L_L         = OBJ.car_tab{"L_L"        , "Value"};
end
%--------------------------------------------------------------------------
function OBJ = LOC_assign_ICs(OBJ)

% NOTE:  the order of teh DOF vector is:
%          1      2    3     4     5     6     7
%      [ y_c, theta, psi, y_bR, y_fR, y_bL, y_fL ]

     OBJ.IC_pos = [ OBJ.IC_tab{"y_c",       "Value"};
                    OBJ.IC_tab{"theta",     "Value"};
                    OBJ.IC_tab{"psi",       "Value"};
                    OBJ.IC_tab{"y_bR",      "Value"};
                    OBJ.IC_tab{"y_fR",      "Value"};
                    OBJ.IC_tab{"y_bL",      "Value"};
                    OBJ.IC_tab{"y_fL",      "Value"};  ];
   
     OBJ.IC_vel = [ OBJ.IC_tab{"y_dot_c",   "Value"};
                    OBJ.IC_tab{"theta_dot", "Value"};
                    OBJ.IC_tab{"psi_dot",   "Value"};
                    OBJ.IC_tab{"y_dot_bR",  "Value"};
                    OBJ.IC_tab{"y_dot_fR",  "Value"};
                    OBJ.IC_tab{"y_dot_bL",  "Value"};
                    OBJ.IC_tab{"y_dot_fL",  "Value"};  ];
     
     
end
%--------------------------------------------------------------------------

