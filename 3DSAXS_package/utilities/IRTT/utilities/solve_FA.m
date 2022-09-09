function FA=solve_FA(tensor,which_eig,thres)
    if nargin<3
        thres=0;
        if nargin<2
            which_eig=1;
        end
    end
    which_eig=4-which_eig;
    nx=size(tensor,1);
    ny=size(tensor,2);
    nz=size(tensor,3);
    FA=zeros(nx,ny,nz);
    for i=1:nx
        for j=1:ny
            for k=1:nz
                xx=tensor(i,j,k,1); yy=tensor(i,j,k,2); zz=tensor(i,j,k,3); xy=tensor(i,j,k,4); xz=tensor(i,j,k,5); yz=tensor(i,j,k,6);
                T=[xx,xy,xz;xy,yy,yz;xz,yz,zz];
                [eigvecs,eigvals]=eig(T);
                eigvals=abs(eigvals);
                if eigvals(which_eig,which_eig)>thres
                    FA(i,j,k)=sqrt(1/2*((eigvals(1,1)-eigvals(2,2))^2+(eigvals(2,2)-eigvals(3,3))^2+(eigvals(3,3)-eigvals(1,1))^2)/(eigvals(1,1)^2+eigvals(2,2)^2+eigvals(3,3)^2));
                end
            end
        end
    end
end
