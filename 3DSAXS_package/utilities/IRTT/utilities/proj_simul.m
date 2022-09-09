function projection=proj_simul(tomogram,R,projsize)
%     load('tomo_myelin.mat','tomotrans');
    tomo = tomogram;
    modelsize = size(tomo);
    if nargin < 3
        projsize = [modelsize(1), modelsize(2)];
    end
    xgrid = (-modelsize(1)/2+1 ):(modelsize(1)/2);
    ygrid = (-modelsize(2)/2+1 ):(modelsize(2)/2);
    zgrid =  (-modelsize(3)/2+1) : (modelsize(3)/2);
    [Y,X,Z] = meshgrid(ygrid, xgrid, zgrid);
    
    Xp = R(2,2)*X + R(2,1)*Y + R(2,3)*Z + projsize(1)/2;
    Yp = R(1,2)*X + R(1,1)*Y + R(1,3)*Z + projsize(2)/2;
%     Xp = R(2,2)*X + R(2,1)*Y + R(1,3)*Z + projsize(1)/2;
%     Yp = R(2,1)*X + R(1,1)*Y + R(2,3)*Z + projsize(2)/2;

    projection=projection_sim(tomo,modelsize,projsize,Xp,Yp); 
end