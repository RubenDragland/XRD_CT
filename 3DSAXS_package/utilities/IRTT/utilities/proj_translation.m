function out=proj_translation(proj,delta)
% Translate and interpolate a projection according to delta=[dx,dy].
    proj=squeeze(proj);
    sizes=size(proj);
    out=zeros(sizes);
    for x=1:sizes(1)
        for y=1:sizes(2)
            xnew=x+delta(1);
            ynew=y+delta(2);
            x2=fix(xnew);
            y2=fix(ynew);
            if (x2>0)&&(x2<=sizes(1))&&(y2>0)&&(y2<=sizes(2))
                out(x2,y2)=out(x2,y2)+proj(x,y)*(x2+1-xnew)*(y2+1-ynew);
            end
            y2=y2+1;
            if (x2>0)&&(x2<=sizes(1))&&(y2>0)&&(y2<=sizes(2))
                out(x2,y2)=out(x2,y2)+proj(x,y)*(x2+1-xnew)*(ynew-y2+1);
            end
            x2=x2+1;
            if (x2>0)&&(x2<=sizes(1))&&(y2>0)&&(y2<=sizes(2))
                out(x2,y2)=out(x2,y2)+proj(x,y)*(xnew-x2+1)*(ynew-y2+1);
            end
            y2=y2-1;
            if (x2>0)&&(x2<=sizes(1))&&(y2>0)&&(y2<=sizes(2))
                out(x2,y2)=out(x2,y2)+proj(x,y)*(xnew-x2+1)*(y2+1-ynew);
            end
        end
    end
end