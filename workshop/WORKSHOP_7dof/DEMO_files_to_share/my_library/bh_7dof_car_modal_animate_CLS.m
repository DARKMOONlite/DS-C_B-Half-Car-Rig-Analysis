classdef bh_7dof_car_modal_animate_CLS
    properties(SetAccess = protected)
        hf
        ha_AB_XY_RIGHT
        ha_DC_XY_LEFT
        ha_BC_YZ_FRONT
        ha_AD_YZ_BACK
    end
%==========================================================================    
    properties
        pause_secs = 0;
    end
%==========================================================================    
   properties(SetAccess = protected)
       view_obj_vec (:,4)  = bh_1view_animate_CLS.empty(0,4)     
   end
    
%==========================================================================    
    methods
        
        %------------------------------------------------------------------
        % ALLOWED calling syntax
        %  >>  func()
        %  >>  func(hf)
        function OBJ = bh_7dof_car_modal_animate_CLS(varargin)
            switch(nargin)
                case 0
                    OBJ.hf               = figure;
                    OBJ.ha_AB_XY_RIGHT   = subplot(2,2,1);
                    OBJ.ha_DC_XY_LEFT    = subplot(2,2,3);
                    OBJ.ha_BC_YZ_FRONT   = subplot(2,2,2);
                    OBJ.ha_AD_YZ_BACK    = subplot(2,2,4);
                case 1
                    OBJ.hf               = varargin{1};
                    OBJ.ha_AB_XY_RIGHT   = subplot(2,2,1);
                    OBJ.ha_DC_XY_LEFT    = subplot(2,2,3);
                    OBJ.ha_BC_YZ_FRONT   = subplot(2,2,2);
                    OBJ.ha_AD_YZ_BACK    = subplot(2,2,4);
                otherwise
                    error("###_ERROR:  UNknown usage mode");
            end
            
            % create the View Objects
            OBJ.view_obj_vec(1) = bh_1view_animate_CLS("side");
            OBJ.view_obj_vec(1).Color_w1  = "blue";
            OBJ.view_obj_vec(1).Color_w2  = "red";
            OBJ.view_obj_vec(1).Name_view = "AB_XY_RIGHT";
            OBJ.view_obj_vec(1) =  OBJ.view_obj_vec(1).birth();

            OBJ.view_obj_vec(2) = bh_1view_animate_CLS("front");
            OBJ.view_obj_vec(2).Color_w1  = "red";
            OBJ.view_obj_vec(2).Color_w2  = [195, 177, 225]/255;  % purple pastel
            OBJ.view_obj_vec(2).Name_view = "BC_YZ_FRONT";
            OBJ.view_obj_vec(2) =  OBJ.view_obj_vec(2).birth();

            OBJ.view_obj_vec(3) = bh_1view_animate_CLS("side");
            OBJ.view_obj_vec(3).Color_w1  = [169, 169, 169]/255;  % grey
            OBJ.view_obj_vec(3).Color_w2  = [195, 177, 225]/255;  % purple pastel
            OBJ.view_obj_vec(3).Name_view = "DC_XY_LEFT";
            OBJ.view_obj_vec(3) = OBJ.view_obj_vec(3).birth();

            OBJ.view_obj_vec(4) = bh_1view_animate_CLS("front");
            OBJ.view_obj_vec(4).Color_w1  = "blue";
            OBJ.view_obj_vec(4).Color_w2  = [169, 169, 169]/255;  % grey
            OBJ.view_obj_vec(4).Name_view = "AD_YZ_BACK";
            OBJ.view_obj_vec(4) = OBJ.view_obj_vec(4).birth();           
        end
        %------------------------------------------------------------------
        function animate(OBJ, eVec_table, mode_ID)
           arguments
             OBJ
             eVec_table (7,7) table
             mode_ID    (1,1) double {mustBeInRange(mode_ID, 1, 7)}             
           end
                     
           % extract the specific dof for the 4 views
           dof_evec(:,1) = eVec_table{["y_c","theta","y_bR","y_fR"], mode_ID}; % "AB_XY_RIGHT"
           dof_evec(:,2) = eVec_table{["y_c","psi",  "y_fR","y_fL"], mode_ID}; % "BC_YZ_FRONT"
           dof_evec(:,3) = eVec_table{["y_c","theta","y_bL","y_fL"], mode_ID}; % "DC_XY_LEFT"
           dof_evec(:,4) = eVec_table{["y_c","psi",  "y_bR","y_bL"], mode_ID}; % "AD_YZ_BACK"
           
           % create a simple sine wave to drive the animation
           y_list       = sin(2*pi*1*[0:0.05:2]);
           ang_deg_list = 15*y_list;
           ang_rad_list = deg2rad(ang_deg_list);
           for kk=1:length(y_list)
               
             x1_dof = dof_evec(:,1) .* [y_list(kk);ang_rad_list(kk);y_list(kk);y_list(kk);];
             x2_dof = dof_evec(:,2) .* [y_list(kk);ang_rad_list(kk);y_list(kk);y_list(kk);];     
             x3_dof = dof_evec(:,3) .* [y_list(kk);ang_rad_list(kk);y_list(kk);y_list(kk);];     
             x4_dof = dof_evec(:,4) .* [y_list(kk);ang_rad_list(kk);y_list(kk);y_list(kk);];                    
               
             OBJ.view_obj_vec(1) = OBJ.view_obj_vec(1).update_poly(x1_dof);
             OBJ.view_obj_vec(2) = OBJ.view_obj_vec(2).update_poly(x2_dof);     
             OBJ.view_obj_vec(3) = OBJ.view_obj_vec(3).update_poly(x3_dof);     
             OBJ.view_obj_vec(4) = OBJ.view_obj_vec(4).update_poly(x4_dof);         

             %axes(hax)
             OBJ.view_obj_vec(1).plot(OBJ.ha_AB_XY_RIGHT); 
             OBJ.view_obj_vec(2).plot(OBJ.ha_BC_YZ_FRONT);      
             OBJ.view_obj_vec(3).plot(OBJ.ha_DC_XY_LEFT);
             OBJ.view_obj_vec(4).plot(OBJ.ha_AD_YZ_BACK);     %fprintf("\n ... %3d", kk)
             drawnow             
             
           end %  for kk=1:length(y_list)
           
        end
        %------------------------------------------------------------------
    end
end
%_#########################################################################
%_ LOCAL SUBFUNCTIONS
%_#########################################################################

