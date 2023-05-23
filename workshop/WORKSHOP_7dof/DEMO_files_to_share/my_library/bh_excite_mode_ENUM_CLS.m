classdef bh_excite_mode_ENUM_CLS < Simulink.IntEnumType
    enumeration
        MODAL_Y_BOUNCE (1)
        MODAL_ROLL     (2)
        MODAL_PITCH    (3)
        MODAL_FRONT_OUT_PHASE (4)
        MODAL_FRONT_IN_PHASE (5)
        MODAL_BACK_OUT_PHASE (6)
        MODAL_BACK_IN_PHASE  (7)
        %-----------------------------
        SINE   (8)
        SINE_ROLL     (9)
        SINE_PITCH    (10)
        SINE_FRONT_OUT_PHASE (11)
        SINE_FRONT_IN_PHASE (12)
        SINE_BACK_OUT_PHASE (13)
        SINE_BACK_IN_PHASE  (14)
        
    end
 %=========================================================================   
    methods
        function tf_out = isSineMember(OBJ)
            tmp_sine_list = [ ...
                    "SINE";
                    "SINE_ROLL";     
                    "SINE_PITCH";    
                    "SINE_FRONT_OUT_PHASE"; 
                    "SINE_FRONT_IN_PHASE"; 
                    "SINE_BACK_OUT_PHASE"; 
                    "SINE_BACK_IN_PHASE";  ];
             tmp_L = ismember(string(OBJ), tmp_sine_list);
             
             tf_out = any(tmp_L);
        end
        
    end
end
