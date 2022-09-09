%CONVERT_TRANSFORM translate between orientation definitions
% trfData = convert_transform(orientation, typeSource, typeTarget)
% These transformations are applied from right to left, e.g. in case of
% 'transpose_fliplr_flipud' and an orientation vector [1 0 1], the vector 
% represents the transform data_mod = transpose(flipud(data)).
% 
% ** orientation            orientation vector for input type typeSource
% ** typeSource             orientation type; 'transpose_fliplr_flipud','transpose_rot90' or 'fliplr_rot90'
% ** typeTarget             orientation type; 'transpose_fliplr_flipud','transpose_rot90' or 'fliplr_rot90'
%
% Examples:
%       trfData = math.convert_transform([1 0], 'fliplr_rot90','transpose_fliplr_flipud');
%       trfData = math.convert_transform([0 1 0], 'transpose_fliplr_flipud','fliplr_rot90');
% 
% see also: math.applyTransform

% Academic License Agreement
%
% Source Code
%
% Introduction 
% •	This license agreement sets forth the terms and conditions under which the PAUL SCHERRER INSTITUT (PSI), CH-5232 Villigen-PSI, Switzerland (hereafter "LICENSOR") 
%   will grant you (hereafter "LICENSEE") a royalty-free, non-exclusive license for academic, non-commercial purposes only (hereafter "LICENSE") to use the cSAXS 
%   ptychography MATLAB package computer software program and associated documentation furnished hereunder (hereafter "PROGRAM").
%
% Terms and Conditions of the LICENSE
% 1.	LICENSOR grants to LICENSEE a royalty-free, non-exclusive license to use the PROGRAM for academic, non-commercial purposes, upon the terms and conditions 
%       hereinafter set out and until termination of this license as set forth below.
% 2.	LICENSEE acknowledges that the PROGRAM is a research tool still in the development stage. The PROGRAM is provided without any related services, improvements 
%       or warranties from LICENSOR and that the LICENSE is entered into in order to enable others to utilize the PROGRAM in their academic activities. It is the 
%       LICENSEE’s responsibility to ensure its proper use and the correctness of the results.”
% 3.	THE PROGRAM IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR 
%       A PARTICULAR PURPOSE AND NONINFRINGEMENT OF ANY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS. IN NO EVENT SHALL THE LICENSOR, THE AUTHORS OR THE COPYRIGHT 
%       HOLDERS BE LIABLE FOR ANY CLAIM, DIRECT, INDIRECT OR CONSEQUENTIAL DAMAGES OR OTHER LIABILITY ARISING FROM, OUT OF OR IN CONNECTION WITH THE PROGRAM OR THE USE 
%       OF THE PROGRAM OR OTHER DEALINGS IN THE PROGRAM.
% 4.	LICENSEE agrees that it will use the PROGRAM and any modifications, improvements, or derivatives of PROGRAM that LICENSEE may create (collectively, 
%       "IMPROVEMENTS") solely for academic, non-commercial purposes and that any copy of PROGRAM or derivatives thereof shall be distributed only under the same 
%       license as PROGRAM. The terms "academic, non-commercial", as used in this Agreement, mean academic or other scholarly research which (a) is not undertaken for 
%       profit, or (b) is not intended to produce works, services, or data for commercial use, or (c) is neither conducted, nor funded, by a person or an entity engaged 
%       in the commercial use, application or exploitation of works similar to the PROGRAM.
% 5.	LICENSEE agrees that it shall make the following acknowledgement in any publication resulting from the use of the PROGRAM or any translation of the code into 
%       another computing language:
%       "Data processing was carried out using the cSAXS ptychography MATLAB package developed by the Science IT and the coherent X-ray scattering (CXS) groups, Paul 
%       Scherrer Institut, Switzerland."
%
% Additionally, any publication using the package, or any translation of the code into another computing language should cite for difference map:
% P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, High-resolution scanning X-ray diffraction microscopy, Science 321, 379–382 (2008). 
%   (doi: 10.1126/science.1158573),
% for maximum likelihood:
% P. Thibault and M. Guizar-Sicairos, Maximum-likelihood refinement for coherent diffractive imaging, New J. Phys. 14, 063004 (2012). 
%   (doi: 10.1088/1367-2630/14/6/063004),
% for mixed coherent modes:
% P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
% and/or for multislice:
% E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089–29108 (2016). 
%   (doi: 10.1364/OE.24.029089).
% 6.	Except for the above-mentioned acknowledgment, LICENSEE shall not use the PROGRAM title or the names or logos of LICENSOR, nor any adaptation thereof, nor the 
%       names of any of its employees or laboratories, in any advertising, promotional or sales material without prior written consent obtained from LICENSOR in each case.
% 7.	Ownership of all rights, including copyright in the PROGRAM and in any material associated therewith, shall at all times remain with LICENSOR, and LICENSEE 
%       agrees to preserve same. LICENSEE agrees not to use any portion of the PROGRAM or of any IMPROVEMENTS in any machine-readable form outside the PROGRAM, nor to 
%       make any copies except for its internal use, without prior written consent of LICENSOR. LICENSEE agrees to place the following copyright notice on any such copies: 
%       © All rights reserved. PAUL SCHERRER INSTITUT, Switzerland, Laboratory for Macromolecules and Bioimaging, 2017. 
% 8.	The LICENSE shall not be construed to confer any rights upon LICENSEE by implication or otherwise except as specifically set forth herein.
% 9.	DISCLAIMER: LICENSEE shall be aware that Phase Focus Limited of Sheffield, UK has an international portfolio of patents and pending applications which relate 
%       to ptychography and that the PROGRAM may be capable of being used in circumstances which may fall within the claims of one or more of the Phase Focus patents, 
%       in particular of patent with international application number PCT/GB2005/001464. The LICENSOR explicitly declares not to indemnify the users of the software 
%       in case Phase Focus or any other third party will open a legal action against the LICENSEE due to the use of the program.
% 10.	This Agreement shall be governed by the material laws of Switzerland and any dispute arising out of this Agreement or use of the PROGRAM shall be brought before 
%       the courts of Zürich, Switzerland. 

function trfData = convert_transform(orientation, typeSource, typeTarget)

if strcmpi(typeSource, typeTarget)==1
    trfData = orientation;
    return
end

switch lower(typeSource)
    case 'transpose_rot90'
        if length(orientation)~=2
            error('''%s'' expects an orientation vector of length %u.', typeSource, 2);
        end
        orientation(1) = mod(orientation(1), 2);
        orientation(2) = mod(orientation(2), 4);
        
        switch lower(typeTarget)
            case 'fliplr_rot90'
                trfData = [0 orientation(2)];
                if all(orientation==0)
                    return
                end
                if orientation(1)==1
                    trfData(1) = trfData(1)+1;
                    trfData(2) = orientation(2)+3;
                end
                trfData(1) = mod(trfData(1), 2);
                trfData(2) = mod(trfData(2), 4);
                
            case 'transpose_fliplr_flipud'
                trfData = [0 0 0];
                if all(orientation==0)
                    return
                end
                if orientation(2)==1
                    trfData(1:2) = trfData(1:2) + 1;
                elseif orientation(2)==2
                     trfData(2:3) = trfData(2:3) + 1;
                elseif orientation(2)==3
                    trfData(3) = trfData(3) + 1;
                    trfData(1) = trfData(1) + 1;
                end
                if orientation(1)==1
                    trfData(1) = trfData(1) + 1;
                end
                trfData = mod(trfData, 2);
                
        end
        
        
    case 'fliplr_rot90'
        if length(orientation)~=2
            error('''%s'' expects an orientation vector of length %u.', typeSource, 2);
        end
        orientation(1) = mod(orientation(1), 2);
        orientation(2) = mod(orientation(2), 4);
        
        switch lower(typeTarget)
            case 'transpose_rot90'
                trfData = [0 orientation(2)];
                if all(orientation==0)
                    return
                end
                if orientation(1)==1
                    trfData = trfData + 1;
                end
                trfData(1) = mod(trfData(1), 2);
                trfData(2) = mod(trfData(2), 4);
            case 'transpose_fliplr_flipud'
                trfData = [0 0 0];
                if all(orientation==0)
                    return
                end
                if orientation(1)==1
                    trfData(2) = 1;
                end
                if orientation(2)==1
                    trfData(1:2) = trfData(1:2) + 1;
                elseif orientation(2)==2
                     trfData(2:3) = trfData(2:3) + 1;
                elseif orientation(2)==3
                    trfData(3) = trfData(3) + 1;
                    trfData(1) = trfData(1) + 1;
                end
                if orientation(1)==1&&mod(orientation(2),2)==1
                    trfData(2:3) = trfData(2:3)-1;
                end
                trfData(1) = mod(trfData(1), 2);
                trfData(2) = mod(trfData(2), 2);
                trfData(3) = mod(trfData(3), 2);
                
        end
        


    case 'transpose_fliplr_flipud'
        if length(orientation)~=3
            error('''%s'' expects an orientation vector of length %u.', typeSource, 3);
        end
        switch lower(typeTarget)
            
            case 'transpose_rot90'
                trfData = [0 0];
                if all(orientation==0)
                    return
                end
                if orientation(1) == 0                    
                    trfData = trfData+1-prod(orientation(2:3));
                    trfData(2) = trfData(2)+2*orientation(3);
                    
                elseif orientation(1) == 1
                    trfData(1) = trfData(1) + (orientation(2)==orientation(3));
                    trfData(2) = trfData(2) + 3*orientation(3) + orientation(2)-2*prod(orientation(2:3));
                end

                % reduce number of rot90
                trfData(2) = mod(trfData(2),4);
                
                % reduce number of permutations
                trfData(1) = mod(trfData(1),2);

 
            case 'fliplr_rot90'
                trfData = [0 0];
                if all(orientation==0)
                    return
                end
                if orientation(1) == 0               
                    trfData = trfData+1-prod(orientation(2:3));
                    trfData(2) = trfData(2)+2*orientation(3);
                    
                elseif orientation(1) == 1
                    trfData(1) = trfData(1) + (orientation(2)==orientation(3));
                    trfData(2) = trfData(2) + 3*orientation(3) + orientation(2)-2*prod(orientation(2:3));
                end

                % reduce number of rot90
                trfData(2) = mod(trfData(2),4);
                
                % reduce number of permutations
                trfData(1) = mod(trfData(1),2);

                if trfData(1)==1
                    trfData(2) = trfData(2) + 3;
                end
                trfData(1) = mod(trfData(1), 2);
                trfData(2) = mod(trfData(2), 4);
        end
        
        
end
end




