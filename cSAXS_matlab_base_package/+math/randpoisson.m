% function out=randpoisson(inarray,thresh);
% outputs an array of poisson-distributed numbers with mean equal to inarray
% For inarray values above the threshold thresh (default=32),
%  use a quick-and-dirty version of the gaussian method,
%  but with negatives clipped to zero
% J.R. Fienup 10/22/99

function out=randpoisson(inarray,thresh)

if nargin <  1
    error('Requires at least one input argument.');
end

if exist('thresh', 'var')~=1, thresh=32; end

out=inarray;
% High-count pixels - use Gaussian approach
gtthresh=find(inarray>thresh);
if ~isempty(gtthresh)
    out(gtthresh)=inarray(gtthresh) + sqrt(inarray(gtthresh)).*randn(size(inarray(gtthresh)));
    out(gtthresh)=round(max(0,out(gtthresh)));
end
% Low-count pixels - this goes into the counting-experiment method

ltthresh=find(inarray<=thresh);
if ~isempty(ltthresh)
    lamda=inarray(ltthresh); % segregate low-value pixels to speed computation
    % Now dealing with a 1-D column vector that will merge into n-D array out later on
    %Initialize r to zero.
    r = zeros(size(lamda));  % output array for ltthresh pixels
    p = zeros(size(lamda));
    ind = true(size(lamda));
    
    while any(ind)
        p(ind) = p(ind) - log(rand(length(ind),1));  % note, do repeatedly calculate over all of lamda
        ind = find(p < lamda); % Q: does this k index over
        r(ind) = r(ind) + 1;
    end
    
    
    % Return NaN if lamda is not positive -- to do this, un-comment what follows (gives zero now).
    
    % 	tmp = NaN;
    % 	if any(any(any(lamda <= 0)));
    % 	    if prod(size(lamda) == 1),   % i.e., a single pixel?
    % 	        r = tmp(ones(size(lamda)));
    % 	    else
    % 	        k = find(lamda <= 0);
    % 	        r(k) = tmp(ones(size(k)));
    % 	    end
    % 	end
    
    out(ltthresh)=r;  % Merge low-value-pixel results with large-value-pixel results
    
end

end

