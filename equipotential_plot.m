function [xx, yy, pot_height] = equipotential_plot(varargin)

   for i_=1:2:length(varargin)
        switch varargin{i_}
            case 'V'
                V = varargin{i_+1};
            case 'V0'
                V0 = varargin{i_+1};
            case 'z'
                z = varargin{i_+1};
            case 'X'
                X = varargin{i_+1};
            case 'Y'
                Y = varargin{i_+1};
            otherwise
                warning(['Unrecognised input ' num2str((i_+1)/2)])
        end
    end

    if ~exist('V0', 'var')
        V0 = 0;
    end
    %disp("V0 = " + num2str(V0,'%.2f'))

    n1 = size(V, 1);
    n2 = size(V, 2);

    % Create an equipotential plot for V=V0, use linear interpolation to get
    % closer to the z value of the V=V0 point.
    inds0 = zeros(n1, n2);
    for ia1=1:n1
        for ia2=1:n2
            %for indz=1:length(z)-1
            %    if V(ia1,ia2,indz) > V0 && V(ia1,ia2,indz+1) < V0
            for indz=length(z)-1:-1:2
                if V(ia1,ia2,indz) < V0 && V(ia1,ia2,indz-1) > V0
                %if V(ia1,ia2,indz+1) < V0
                    m = (V(ia1,ia2,indz) - V(ia1, ia2, indz+1))/(indz - (indz+1));
                    c = V(ia1, ia2,indz) - m*indz;
                    inds0(ia1,ia2) = -c/m;
                    %inds0(ia1,ia2) = indz;
                    break;
                end
            end
        end
    end
    pot_height = z(1) + inds0*(z(2)-z(1));

    %figure
    %f = figure('visible','off');
    potMin = min(pot_height,[],"all");
    potMax = max(pot_height,[],"all");
    %disp("Min = " + num2str(potMin) + ", max = " + num2str(potMax));
    surf(X, Y, pot_height)
    zMin = min(z,[],"all");
    zMax = max(z,[],"all");
    zlim([zMin,zMax])
    xlabel('x/Å')
    ylabel('y/Å')
    title(['Equipotential V=' num2str(V0,'%.2f') ' meV, used in simulation'])
    %fontsize(gcf,scale=2)
    %daspect([1 1 1])
    %colormap(plasma)
    if(potMin - potMax > 1e-6)
      clim([potMin potMax])
    end
end
