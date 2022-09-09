%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
%*-----------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computing language this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language the authors and institution should be acknowledged 
%   in written form in the publication: “Data processing was carried out 
%   using the “cSAXS software package” developed by the CXS group,
%   Paul Scherrer Institut, Switzerland.” 
%   Variations on the latter text can be incorporated upon discussion with 
%   the CXS group if needed to more specifically reflect the use of the package 
%   for the published work.
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.



lfile = './license.txt';                        % txt file containing the new license
package_path = '/das/work/p16/p16602/code_testing_env3/cSAXS_matlab_ptycho/'; % path to the git repository
start_string = 'Academic License Agreement';    % (sub)string indicating the start of the license
end_string = 'courts of Zürich, Switzerland';   % (sub)string indicating the end of the license



[~, out] = system(sprintf('find %s -name "*.m"',package_path));

fnameList = strsplit(strtrim(out));

for fileID=1:numel(fnameList)
    if exist(fnameList{fileID}, 'file')
        fname = fnameList{fileID};
        fnameTmp = [fname '.tmp'];
    else
        continue;
    end
    
    fid = fopen(fname, 'r', 'n', 'UTF-8');
    fout = fopen([fname '.tmp'], 'w', 'n', 'UTF-8');
    
    content = {};
    licenseStart = 0;
    tline = fgetl(fid);
    while ischar(tline)
        if contains(tline, start_string)
            licenseStart = 1;
        elseif contains(tline, end_string)
            licenseStart = 0;
            printLicense(fout, lfile);
        elseif ~licenseStart
            fprintf(fout, '%s\n', tline);
        end
        tline = fgetl(fid);
    end
    
    movefile(fnameTmp, fname);
    
end


function printLicense(fout, lfile)
flic = fopen(lfile, 'r', 'n', 'UTF-8');
tline = fgetl(flic);
while ischar(tline)
    fprintf(fout, '%s\n', tline);
    tline = fgetl(flic);
end

end