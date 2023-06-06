function [data,init_var] = interpretdata(LVM_DATA)
    data.rawdata = LVM_DATA.Segment1.data
    data.catagoriseddata.time = LVM_DATA.Segment1.data(:,1)
    data.catagoriseddata.LVDT1 = LVM_DATA.Segment1.data(:,2)
    data.catagoriseddata.LVDT2 = LVM_DATA.Segment1.data(:,3)
    data.catagoriseddata.LVDT3 = LVM_DATA.Segment1.data(:,4)
    data.catagoriseddata.LVDT4 = LVM_DATA.Segment1.data(:,5)
    
    
    [~,peaks] = findpeaks(LVM_DATA.Segment1.data(:,2))
    
    init_var = LVM_DATA.Segment1.data(peaks(1),2:5)
end