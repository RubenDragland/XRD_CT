function [theta_out,phi_out,a_out] = Tensor_to_SH(tomotensor)
    sampsize=size(tomotensor);
    phi_out=zeros(sampsize(1:3));
    theta_out=zeros(sampsize(1:3));
    a_out=zeros([sampsize(1:3),4]);
    for i=1:sampsize(1)
        for j=1:sampsize(2)
            for k=1:sampsize(3)
                xx = tomotensor(i,j,k,1); 
                yy = tomotensor(i,j,k,2); 
                zz = tomotensor(i,j,k,3); 
                xy = tomotensor(i,j,k,4); 
                xz = tomotensor(i,j,k,5); 
                yz = tomotensor(i,j,k,6);
                T = [xx,xy,xz;xy,yy,yz;xz,yz,zz];
                [eigvecs, eigvals] = eig(T);
                eigvals(eigvals<0)=0;
                if (eigvals(3,3)-eigvals(2,2))>(eigvals(2,2)-eigvals(1,1))
                    eigvals=sqrt(eigvals);
                    a_out(i,j,k,1)=2*sqrt(pi)*(eigvals(1,1)+eigvals(2,2)+eigvals(3,3))/3;
                    a_out(i,j,k,2)=2/3*sqrt(pi/5)*(-eigvals(1,1)-eigvals(2,2)+2*eigvals(3,3));
                    which_eig=3;
                else
                    eigvals=sqrt(eigvals);
                    a_out(i,j,k,1)=2*sqrt(pi)*(eigvals(1,1)+eigvals(2,2)+eigvals(3,3))/3;
                    a_out(i,j,k,2)=2/3*sqrt(pi/5)*(2*eigvals(1,1)-eigvals(2,2)-eigvals(3,3));
                    which_eig=1;
                end
                v(1:3)=eigvecs([2 1 3],which_eig);
                if v(3)~=0
                    if v(3)>0
                        theta_out(i,j,k)=atan(sqrt(v(1)^2+v(2)^2)/v(3));
                    else
                        theta_out(i,j,k)=atan(sqrt(v(1)^2+v(2)^2)/v(3))+pi;
                    end
                else
                    theta_out(i,j,k)=pi/2;
                end
                if theta_out(i,j,k)~=0
                    phi_out(i,j,k)=cart2pol(v(1)/sin(theta_out(i,j,k)),v(2)/sin(theta_out(i,j,k)));
                end
            end
        end
    end
%     theta_out(:)=3*pi/4;
%     phi_out(:)=pi/2;
    theta_out=reshape(theta_out,1,[]);
    phi_out=reshape(phi_out,1,[]);
    a_out=reshape(a_out,1,[]);
end
                    
                    
        
    
    