function out=proj_from_peak(inputfileformat,firstscan,numlines,peakname,property)
%  proj_from_peak: get a projection from the result of peak fitting
%  fistscan: scan number of the first scan in the projection
%  numlines: number of lines in the projection
%  peakname: name of the peak to use
%  property: =1 to take peak area, =0 to take peak height
    
    if property==1
        p=5;
    else
        p=3;
    end
    
    load(num2str(firstscan,inputfileformat));
    numpoints=eval(['size(' peakname ',1)']);
    out=zeros(numpoints,numlines);
    for j=1:numlines
        load(num2str(firstscan+j-1,inputfileformat));
        peaks=eval(peakname);
        if mod(j,2)==0
            out(:,j)=peaks(:,p);
        else
            out(:,j)=flip(peaks(:,p));
        end
    end
end