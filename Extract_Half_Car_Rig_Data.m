function [averaged_data]=Extract_Half_Car_Rig_Data()

    
    [results] = find_data("LAB TEST FOLDER");
    %% Extracting Relavent Data
    for i =1:size(results,1)
        [data(i),init_var] = interpretdata(lvm_import(fullfile(results(i).folder,results(i).name)));
        folders=regexp(fullfile(results(i).folder,results(i).name),'\','split');
        data_folder_names(i,1:2) = [folders(size(folders,2)-1),folders(size(folders,2)-2)];
        name = regexp(fullfile(results(i).folder,results(i).name),'_','split');
        if size(name,2)>=4
            temp=regexp(name(4),'\.','split');
            data_folder_names(i,3) = cellstr(temp{1}{1});
        else
            data_folder_names(i,3) = {'0'};
        end
        if size(name,2)>=5
            temp=regexp(name(5),'\.','split');
           data_folder_names(i,4)=cellstr(temp{1}{1});
    
        else
            data_folder_names(i,4)={'0'};
        end
        
    end
    for i=1:size(data,2)
    
        data(i).("Test") = data_folder_names(i,1);
        data(i).("Lift_Position") = data_folder_names(i,2);
        data(i).("damping") = data_folder_names(i,3);
        data(i).("mass") = data_folder_names(i,4);
        data(i).("init_var")= init_var;
    end
    %% Sorting Data
    test = struct2table(data);
    test = sortrows(test,["Test","Lift_Position","damping","mass"]);
    data = table2struct(test);
    j=1;
    %% Averaging Data
    for i=1:3:size(test)
        if strcmp(data(i).Test,data(i+1).Test)&&strcmp(data(i).Test,data(i+2).Test)...
            && strcmp(data(i).Lift_Position,data(i+1).Lift_Position) && strcmp(data(i).Lift_Position,data(i+2).Lift_Position)...
            && strcmp(data(i).damping,data(i+1).damping) && strcmp(data(i).damping,data(i+2).damping)...
            && strcmp(data(i).mass,data(i+1).mass) && strcmp(data(i).mass,data(i+2).mass)
        
            averaged_data(j)=data(i);
    
            combined_array = cat(3,data(i).rawdata,data(i+1).rawdata,data(i+2).rawdata);
            mean(combined_array,3);
            averaged_data(j).rawdata=mean(combined_array,3);
    
            averaged_data(j).cdata.LVDT1=averaged_data(j).rawdata(:,2);
            averaged_data(j).cdata.LVDT2=averaged_data(j).rawdata(:,3);
            averaged_data(j).cdata.LVDT3=averaged_data(j).rawdata(:,4);
            averaged_data(j).cdata.LVDT4=averaged_data(j).rawdata(:,5);


            averaged_data(j).cdata.x1 = (averaged_data(j).cdata.LVDT1+averaged_data(j).cdata.LVDT2)/2;
            averaged_data(j).cdata.x2 = averaged_data(j).cdata.LVDT3;
            averaged_data(j).cdata.x3 = averaged_data(j).cdata.LVDT4;
            averaged_data(j).cdata.roll=atan((averaged_data(j).cdata.LVDT1-averaged_data(j).cdata.LVDT2)/1500);
            averaged_data(j).rawdof = [averaged_data(j).cdata.time,averaged_data(j).cdata.x1,averaged_data(j).cdata.roll,averaged_data(j).cdata.x2,averaged_data(j).cdata.x3];
            
            [~,peaks] = findpeaks(averaged_data(j).rawdof(:,2));

            averaged_data(j).init_var = averaged_data(j).rawdof(peaks(1),2:5);

            j=j+1;
        else
            error("data has dissimilar fields")
            pause
        end
    
    end

end



