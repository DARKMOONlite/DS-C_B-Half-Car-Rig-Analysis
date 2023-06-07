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
plot(data(46).catagoriseddata.time,data(46).rawdata(:,2:5))
plot(data(47).catagoriseddata.time,data(47).rawdata(:,2:5))
plot(data(48).catagoriseddata.time,data(48).rawdata(:,2:5))

row = data(46).rawdata(:,2:5);
fftResult = fft(row);
magnitude = abs(fftResult);
plot(magnitude)
plot(fftResult)

Fs = 5000; % Your Sampling frequency
N = 0.3*(row); % Length of your signal
f = (0:N-1)*(Fs/N); % Frequency vector
plot(f, magnitude)
xlabel('Frequency (Hz)')
ylabel('Magnitude')

%%these aren't right
% Assuming these are the additional row indices
rows = [46, 47, 48]; 
% Extract the rows and compute their average
averageRows = mean(data(46).rawdata(rows, 2:5), 1)
