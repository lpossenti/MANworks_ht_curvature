function starting_point=compute_starting_point_position(Y1,Y2,Z, Size_y, Size_z)
%% compute_starting_point_position finds the coordinates of the points from which the pattern networks will start
% 
%               Y1 = Number of points on the odd index planes
%               Y2 = Number of points on the even index planes
%               Z = Half number of planes 
%               Size_y = size of domain along y direction
%               Size_y = size of domain along z direction
% 
% if Y1-Y2=1 the repetion creates an homogeneous grid where the points of a
% plane are situated in the middle of two consecutive points of
% the consecutives planes.
%   es. 5 and 4   ___________
%                | . . . . . |
%                |  . . . .  |
%                | . . . . . |
%                |  . . . .  |
%                |___________|
%
% otherwise the planes are indipendently subdivided
%
%   es. 5 and 3   ___________
%                | . . . . . |
%                |  .  .  .  |
%                | . . . . . |
%                |  .  .  .  |
%                |___________|
%
%
%   Author: Simone Di Gregorio
%   Politecnico di Milano, 09/02/2017
%   Contact: simone.digre@gmail.com  
%%

DeltaZ=Size_z/(Z*2+1);

if abs(Y1-Y2)==1 % FIRST CASE
    
    if Y1 >= Y2
        DeltaY=Size_y/(Y1*2);

        j=1;

        for i=1:(2*Z)
            if mod(i,2)
                for k=1:Y1
                starting_point(j,1)=j;
                starting_point(j,2)=(2*k-1)*DeltaY;
                starting_point(j,3)=i*DeltaZ;
                j=j+1;
                end
            else 
                for k=1:Y2
                starting_point(j,1)=j;
                starting_point(j,2)=(2*k)*DeltaY;
                starting_point(j,3)=i*DeltaZ;
                j=j+1;
                end
            end
        end

    else 
        
        DeltaY=Size_y/(Y2*2);
        j=1;

        for i=1:(2*Z)
            if mod(i,2)
                for k=1:Y1
                starting_point(j,1)=j;
                starting_point(j,2)=(2*k)*DeltaY;
                starting_point(j,3)=i*DeltaZ;
                j=j+1;
                end
            else 
                for k=1:Y2
                starting_point(j,1)=j;
                starting_point(j,2)=(2*k-1)*DeltaY;
                starting_point(j,3)=i*DeltaZ;
                j=j+1;
                end
            end
        end
    
    
    end
    
else % SECOND CASE abs(Y1-Y2)!=1
    
    DeltaY1=Size_y/(Y1*2);
    DeltaY2=Size_y/(Y2*2);
    
    j=1;

    for i=1:(2*Z)
        if mod(i,2)
            for k=1:Y1
            starting_point(j,1)=j;
            starting_point(j,2)=(2*k-1)*DeltaY1;
            starting_point(j,3)=i*DeltaZ;
            j=j+1;
            end
        else 
            for k=1:Y2
            starting_point(j,1)=j;
            starting_point(j,2)=(2*k-1)*DeltaY2;
            starting_point(j,3)=i*DeltaZ;
            j=j+1;
            end
        end
    end

end


