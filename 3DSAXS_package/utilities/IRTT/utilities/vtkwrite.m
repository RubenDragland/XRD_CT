%% vtkwrite(filename,resolution,type,title,data)
% VTKWRITE Writes 3D Matlab array into VTK file format.
%
% Inputs:
%   filename    name and directory of output vtk file
%   resolution  voxel size of the data, can be arbitrary unit
%   type        'SCALARS', 'VECTORS' and 'TENSORS' are supported
%   title       name of the sample, tissue, etc.
%   data        data matrix, with size (nx*ny*nz) for 'SCALARs', (nx*ny*nz*3)
%   for 'VECTORS' and (nx*ny*nz*6) for 'TENSORS'

function vtkwrite(filename,resolution,type,title,data)
% VTKWRITE Writes 3D Matlab array into VTK file format.

%  vtkwrite(filename,resolution,'scalars',title,r) writes a 3D
%  scalar data into VTK file. r is the scalar data matrix.
%
%  vtkwrite(filename,resolution,'vectors',title,v) writes a 3D vector map
%  into VTK file. v is the vector map matrix (nx*ny*nz*3)

    fid = fopen(filename, 'w'); 
    fprintf(fid, '# vtk DataFile Version 2.0\n');
    fprintf(fid, 'VTK from Matlab\n');
    fprintf(fid, 'ASCII\n');
    nx = size(data,1);
    ny = size(data,2);
    nz = size(data,3);
    
    fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
    fprintf(fid, 'DIMENSIONS %d %d %d\n', nx, ny, nz);
    fprintf(fid, 'SPACING %f %f %f\n', resolution, resolution, resolution);
    fprintf(fid, 'ORIGIN 0 0 0\n');
    fprintf(fid, 'POINT_DATA %d\n', nx*ny*nz);
    
    switch upper(type)
        case 'SCALARS' 
            fprintf(fid, ['SCALARS ',title,' float\n']);
            fprintf(fid, 'LOOKUP_TABLE default\n');
            fprintf(fid, '%.2f ', data(:));
        case 'VECTORS'
            fprintf(fid, ['VECTORS ',title,' float\n']);
            output=[reshape(data(:,:,:,1),1,[]);reshape(data(:,:,:,2),1,[]);reshape(data(:,:,:,3),1,[])];
            fprintf(fid, '%.2f ', output(:));
        case 'TENSORS'
            fprintf(fid, ['TENSORS ',title,' float\n']);
            output=[reshape(data(:,:,:,1),1,[]);reshape(data(:,:,:,4),1,[]);reshape(data(:,:,:,5),1,[]);...
                reshape(data(:,:,:,4),1,[]);reshape(data(:,:,:,2),1,[]);reshape(data(:,:,:,6),1,[]);...
                reshape(data(:,:,:,5),1,[]);reshape(data(:,:,:,6),1,[]);reshape(data(:,:,:,3),1,[]);...
                ];
            fprintf(fid, '%.2f ', output(:));
    end
    fclose(fid);
end