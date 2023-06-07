%% Extract Data
data = Extract_Half_Car_Rig_Data()

%% Plot Data
tiledlayout(8,4)
for i=1:size(data,2)
    nexttile
    plot(data(i).catagoriseddata.time,data(i).rawdata(:,2:5))
    title(strcat(data(i).Lift_Position,"\_",data(i).Test,"\_",data(i).damping,"\_",data(i).mass))
end
