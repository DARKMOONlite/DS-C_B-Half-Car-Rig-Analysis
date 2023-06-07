function [data,init_var] = interpretdata(LVM_DATA)
    data.rawdata = LVM_DATA.Segment1.data;
    data.cdata.time = LVM_DATA.Segment1.data(:,1);
    data.cdata.LVDT1 = LVM_DATA.Segment1.data(:,2);
    data.cdata.LVDT2 = LVM_DATA.Segment1.data(:,3);
    data.cdata.LVDT3 = LVM_DATA.Segment1.data(:,4);
    data.cdata.LVDT4 = LVM_DATA.Segment1.data(:,5);


    data.cdata.x1 = (LVM_DATA.Segment1.data(:,2)+LVM_DATA.Segment1.data(:,3))/2;
    data.cdata.x2=data.cdata.LVDT3;
    data.cdata.x3=data.cdata.LVDT4;
    data.cdata.roll=atan((data.cdata.LVDT1-data.cdata.LVDT2)/1500);
    

  
    [~,peaks] = findpeaks(LVM_DATA.Segment1.data(:,2));
    
    init_var = LVM_DATA.Segment1.data(peaks(1),2:5);
end