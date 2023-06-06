clear all
% test_data = lvm_import("test.lvm")
% [extracted_data,init_param]= interpretdata(test_data)
% 
% 
% plot(extracted_data.catagoriseddata.time,extracted_data.rawdata(:,2:5))
% hold on
% yline(init_param)



results = find_data("LAB TEST FOLDER")

for i =1:size(results,1)
    data(i) = interpretdata(lvm_import(fullfile(results(i).folder,results(i).name)))
    folders=regexp(fullfile(results(i).folder,results(i).name),'\','split')
    data_folder_names(i,1:2) = [folders(size(folders,2)-1),folders(size(folders,2)-2)]
    name = regexp(fullfile(results(i).folder,results(i).name),'_','split')
    if size(name,2)>=4
        temp=regexp(name(4),'\.','split')
        data_folder_names(i,3) = cellstr(temp{1}{1})
        
    end
    if size(name,2)>=5
        temp=regexp(name(5),'\.','split')
        data_folder_names(i,4)=cellstr(temp{1}{1})
        
    end

end


plot(data(90).catagoriseddata.time,data(90).rawdata(:,2:5))