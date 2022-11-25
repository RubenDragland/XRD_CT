function [chi,phi] = get_ID15_SASTT_angles()


sample_name = 's3t2';

datadir_int = 'H:/phd_data/id15a_2018_10/int_xrdct_data/fkm/s3t2_rad_2048_az_32/';


flist_proj = dir(sprintf([datadir_int sample_name '*2048*.mat']));
flist_proj = removeLettersFromProjectionName(flist_proj,sample_name);
flist_proj = natsortfiles({flist_proj.name});
% flist_proj = removeInvalidProjectionNames(flist_proj,sample_name);





for ii = 1:length(flist_proj)

    fname = [datadir_int flist_proj{ii}];
    %get rotation angles from file name
    tmp = strsplit(fname,'phi');
    tmp = strsplit(tmp{2},'_');
    phi(ii) = str2double(tmp{1});
    tmp = strsplit(fname,'chi');
    tmp = strsplit(tmp{2},'_');
    chi(ii) = str2double(tmp{1});
end

end

