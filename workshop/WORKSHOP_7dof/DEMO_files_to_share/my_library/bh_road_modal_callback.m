function varargout = bh_road_modal_callback(action_mode, the_gcb)
    arguments
        action_mode (1,1) string
        the_gcb     (1,:) char
    end


    switch(upper(action_mode))
        case "M_EXCITE_CBACK"
               LOC_do_excite_callback(the_gcb);
        case "GET_FHZ_MODAL"
            f_hz_modal = LOC_get_fhz_modal(gcb);
            varargout{1} = f_hz_modal;
        case "EDIT_THIS_FILE"
            edit(mfilename);
        case "POPULATE_TABLE"
            LOC_populate_table(the_gcb);
        otherwise
            error("###_ERROR:  UNknown usage mode.");
    end

end % function bh_road_modal_callback
%_#########################################################################
% LOCAL Subfunctions
%_#########################################################################
function LOC_do_excite_callback(the_gcb)

    tmp    =  Simulink.Mask.get(the_gcb);
    tmp_ex = tmp.getParameter("m_excite_type");
    tmp_f =  tmp.getParameter("f_str");  

    tmp_sine_list = [ ...
            "bh_excite_mode_ENUM_CLS.SINE";
            "bh_excite_mode_ENUM_CLS.SINE_ROLL";     
            "bh_excite_mode_ENUM_CLS.SINE_PITCH";    
            "bh_excite_mode_ENUM_CLS.SINE_FRONT_OUT_PHASE"; 
            "bh_excite_mode_ENUM_CLS.SINE_FRONT_IN_PHASE"; 
            "bh_excite_mode_ENUM_CLS.SINE_BACK_OUT_PHASE"; 
            "bh_excite_mode_ENUM_CLS.SINE_BACK_IN_PHASE";  ];

     tmp_L = ismember(string(tmp_ex.Value), tmp_sine_list);   

     if( any(tmp_L) )
         tmp_f.Enabled = "on";
     else
         tmp_f.Enabled = "off";
     end

end
%==========================================================================
function out_f_hz_modal = LOC_get_fhz_modal(the_gcb)

    tmp    =  Simulink.Mask.get(the_gcb);
    tmp_ex = tmp.getParameter("m_excite_type");

    m_excite_type = eval(tmp_ex.Value);
    m_excite_str  = string(m_excite_type); 
    
    mode_fhz      = [ ...
                          0.9162; % MODE_1
                          1.1836; % MODE_2
                          1.6742; % MODE_3
                          5.6358; % MODE_4
                          5.6435; % MODE_5
                          6.2992; % MODE_6
                          6.3188; % MODE_7
                        ];
                    
    switch(m_excite_str)
    
        case "MODAL_Y_BOUNCE"
            out_f_hz_modal = mode_fhz(1);
        case "MODAL_ROLL"
            out_f_hz_modal = mode_fhz(2);
        case "MODAL_PITCH"
            out_f_hz_modal = mode_fhz(3);
        case "MODAL_FRONT_OUT_PHASE"
            out_f_hz_modal = mode_fhz(4);
        case "MODAL_FRONT_IN_PHASE"
            out_f_hz_modal = mode_fhz(5);
        case "MODAL_BACK_OUT_PHASE"
            out_f_hz_modal = mode_fhz(6);
        case "MODAL_BACK_IN_PHASE"
            out_f_hz_modal = mode_fhz(7);
        otherwise
            f_hz_modal = 0;
    end        
end
%==========================================================================
function  LOC_populate_table(the_gcb)

    mode_fhz      = [ ...
                          0.9162; % MODE_1
                          1.1836; % MODE_2
                          1.6742; % MODE_3
                          5.6358; % MODE_4
                          5.6435; % MODE_5
                          6.2992; % MODE_6
                          6.3188; % MODE_7
                        ];

    mode_str = string(mode_fhz);

    my_rows = join( ...
        [  "{'1', 'BOUNCE',"+"'"+mode_str(1)+"'"+";",  ...
           " '2', 'ROLL',"+"'"+mode_str(2)+"'"+";",  ...   
           " '3', 'PITCH',"+"'"+mode_str(3)+"'"+";",  ...   
           " '4', 'FRONT_OUT_PHASE',"+"'"+mode_str(4)+"'"+";",  ...
           " '5', 'FRONT_IN_PHASE',"+"'"+mode_str(5)+"'"+";",  ...
           " '6', 'BACK_OUT_PHASE',"+"'"+mode_str(6)+"'"+";",  ...
           " '7', 'BACK_IN_PHASE',"+"'"+mode_str(7)+"'"+"}" ]  );

    tmp = Simulink.Mask.get(gcb);
    tmp_T = tmp.getParameter('m_mode_table');   
    tmp_T.Value = my_rows;
   
end