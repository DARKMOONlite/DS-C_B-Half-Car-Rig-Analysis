classdef bh_1view_animate_CLS
    properties(SetAccess = protected)
        View (1,1) string {mustBeMember(View, ["side", "front"])} = "side"
    end
%==========================================================================    
    properties
        Name_view
        Color_car = "green"
        Color_w1  = "blue"
        Color_w2  = "red"
    end
%==========================================================================    
   properties(SetAccess = protected)
       birth_poly_car
       birth_poly_w1
       birth_poly_w2 
       %---------------
       poly_car
       poly_w1
       poly_w2 
   end
    
%==========================================================================    
    methods
        
        function OBJ = bh_1view_animate_CLS(View)
            arguments
              View (1,1) string {mustBeMember(View, ["side", "front"])};
            end
            
            OBJ.View = View;
        end
        %------------------------------------------------------------------
        function OBJ = birth(OBJ)
            
            OBJ.birth_poly_car = birth_poly(OBJ,"car");
            OBJ.birth_poly_w1  = birth_poly(OBJ,"w1");
            OBJ.birth_poly_w2  = birth_poly(OBJ,"w2");   
                 
            OBJ.poly_car       = OBJ.birth_poly_car;
            OBJ.poly_w1        = OBJ.birth_poly_w1;
            OBJ.poly_w2        = OBJ.birth_poly_w2;           
        end
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function OBJ = update_poly(OBJ, x_dof)
           arguments
                OBJ
                x_dof (4,1) double
           end
               
           poly_car = birth_poly(OBJ,"car");
           poly_w1  = birth_poly(OBJ,"w1");
           poly_w2  = birth_poly(OBJ,"w2");   
           
           poly_car = translate(poly_car, [0, x_dof(1)] );
           [cent_x, cent_y]  = centroid(poly_car);
           poly_car = rotate(poly_car, rad2deg(x_dof(2)),  [cent_x, cent_y] );
  
           poly_w1 = translate(poly_w1, [0, x_dof(3)] );
           poly_w2 = translate(poly_w2, [0, x_dof(4)] );    
           
           OBJ.poly_car = poly_car;
           OBJ.poly_w1  = poly_w1;
           OBJ.poly_w2  = poly_w2;
        end
        %------------------------------------------------------------------
        function plot(OBJ, ha)
            
            % NOTE:  the following code produces a 2nd axes
            %    cla(ha)
            %    hold(ha, "on")
            % So I'll us the DELETE function instead
            hold(ha, "on")           
            z = findobj(ha,"Type","Polygon");
            delete(z);
            
            OBJ.poly_car.plot("Parent", ha, "FaceColor", OBJ.Color_car, "FaceAlpha", 1);            
            OBJ.poly_w1.plot("Parent",  ha, "FaceColor", OBJ.Color_w1, "FaceAlpha", 1);
            OBJ.poly_w2.plot("Parent",  ha, "FaceColor", OBJ.Color_w2, "FaceAlpha", 1);
            %hold(ha,"off");
            axis(ha, "equal")
            grid(ha, "on")
            xlim(ha, [-1, 10])
            ylim(ha, [-1,  8])
            title(ha,OBJ.Name_view, "Interpreter","none")
            %drawnow
        end
        %------------------------------------------------------------------        
        %------------------------------------------------------------------
    end
%-*************************************************************************
methods (Access = protected)
     function my_poly = birth_poly(OBJ,Name)
        arguments
            OBJ
            Name (1,1) string {mustBeMember(Name, ["car","w1","w2"])}
        end

         x_squa  = [0, 1, 1, 0]; 
         y_squa  = [0, 0, 1, 1];
         
         x_squa  = [0, 0.5, 0.5, 0]; 
         y_squa  = [0, 0, 1, 1];

         x_circ  = 0.5 + 0.5*cosd(0:1:360);
         y_circ  = 0.5 + 0.5*sind(0:1:360);
                  
        if(OBJ.View=="side")
           x_rect = [0, 7, 7, 0 ];
           y_rect = [0, 0, 1, 1 ];
           x_w    = x_circ;
           y_w    = y_circ;
           
           w2_offset = x_rect(2) - x_rect(1);
        elseif(OBJ.View=="front")
            % front
           x_rect = [0, 5, 5, 0 ];
           y_rect = [0, 0, 1, 1 ];
           
           x_w    = x_squa;
           y_w    = y_squa;
           
           w2_offset = 0.5 + x_rect(2) - x_rect(1);
        else
            error("###_ERROR:  UNknown view specified");
        end

        switch(Name)
            case "w1"
                my_poly            = polyshape(x_w, y_w);
                my_poly            = translate(my_poly, [1, 1] );                    
            case "w2"
                my_poly            = polyshape(x_w, y_w);
                my_poly            = translate(my_poly, [w2_offset, 1] );                    
            case "car"
                my_poly            = polyshape(x_rect, y_rect);
                my_poly            = translate(my_poly, [1, 5] );                    
        end
    end
    %------------------------------------------------------------------
    %------------------------------------------------------------------  
end % methods (Access = protected)
%-*************************************************************************    
    
end
%_#########################################################################
%_ LOCAL SUBFUNCTIONS
%_#########################################################################

