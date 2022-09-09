% SENDSMS send message to given phone number(s)
% ** recipient          phone number(s);
% ** msg                message
% 
% *optional*
% ** binary_path        path to sendmail binary file
% ** server             server address
% ** sleep              sleep time between messages
% ** logfile            logfile+path to keep track of past activities
%
% EXAMPLES
%       sendSMS('0041123456789', 'done!');
%       sendSMS({'076123456', '076987654'}, 'done!')
%       sendSMS({'076123456', '076987654'}, {'message1', 'message2'})
%       sendSMS({'076123456', '076987654'}, 'done!', 'sleep', 3600)
%
%

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
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

function sendSMS(recipient, msg, varargin)

persistent last_call

par = inputParser;
par.addParameter('binary_path', '/usr/sbin/sendmail', @ischar) 
par.addParameter('server', '@sms.switch.ch', @ischar)
par.addParameter('sleep', 100, @isnumeric)
par.addParameter('logfile', [], @ischar)
par.parse(varargin{:})

var = par.Results;



if isempty(var.logfile)
    if isempty(last_call)
        last_call = datetime('now');
        
    elseif datetime('now')-last_call < seconds(var.sleep)
        return
        
    else
        last_call = datetime('now');
    end
else
    if exist(var.logfile, 'file')
        fid = fopen(var.logfile);
        while ~feof(fid)
            tline = fgetl(fid);
        end
        fclose(fid);
        tline = datetime(tline);
        if datetime('now')-tline < seconds(var.sleep)
            return
        else
            system(['echo "' datestr(datetime('now')) '" >> ' var.logfile]);
        end
    else
        system(['echo "' datestr(datetime('now')) '" >> ' var.logfile]);
    end
end
    



%% input checks
if ~iscell(recipient)
    tmp = recipient;
    recipient = [];
    recipient{1} = tmp;
    clear tmp;
end

if iscell(msg)
    multi_msg = true;
else
    multi_msg = false;
end


%% compile and send message
if ~multi_msg
    for ii=1:numel(recipient)
        if ~ischar(recipient{ii})
            error('The phone number has to be given as a string.')
        end
        msg_full = sprintf('%s\n\n%s', ['To:' recipient{ii} var.server], msg);
        system(['echo -e "' msg_full '" | ' var.binary_path ' -t']);
        
    end
else
    assert(numel(recipient)==numel(msg), 'Number of recipients and messages does not match!');
    for ii=1:numel(recipient)
        if ~ischar(recipient{ii})
            error('The phone number has to be given as a string.')
        end
        msg_full = sprintf('%s\n\n%s', ['To:' recipient{ii} var.server], msg{ii});
        system(['echo -e "' msg_full '" | ' var.binary_path ' -t']);
        
    end
end



end

