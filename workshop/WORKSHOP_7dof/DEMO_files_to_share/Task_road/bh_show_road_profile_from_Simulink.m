function bh_show_road_profile_from_Simulink(the_BDROOT)

    speed_mps = str2num(get_param(the_BDROOT+"/Car speed (m//sec)","Value"));
    Hy          = str2num(get_param(the_BDROOT+"/Hy: Hump Height (m)","Value"));
    Lx          = str2num(get_param(the_BDROOT+"/Lx: Hump Width (m)","Value"));
    x_offset    = str2num(get_param(the_BDROOT+"/x_offset: Amount of Flat road before FIRST hump (m)","Value"));
    NUM_PERIODS = str2num(get_param(the_BDROOT+"/Number of HUMPS","Value"));
    P           = str2num(get_param(the_BDROOT+"/P: Distance between Humps (m)","Value"));

    % create SYS object
    OBJ = bh_road_hump_A_SYSO("Hy",Hy, "Lx",Lx, "P",P, ...
                              "x_offset",x_offset, "NUM_PERIODS",NUM_PERIODS);
    % plot stuff                                       
    figure;
    hax(1) = subplot(2,1,1);
             OBJ.plot_X_domain("Parent", hax(1));
    %  
    hax(2) = subplot(2,1,2);
    OBJ.plot_t_domain(speed_mps, "units", "m/s", "Parent", hax(2) );

end