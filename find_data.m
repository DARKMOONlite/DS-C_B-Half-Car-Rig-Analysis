function [file_list] = find_data(folder_name)

file_list = dir(fullfile(pwd,folder_name,"**\*.lvm"))
file_list = file_list(~[file_list.isdir]);

end