function bar_werror(x,y,error_margin,varargin)
    %inputs x is x data or labels along xaxis
    %       y is y data
    %       error_margin is the error margin for the y data
    %       varargin is the optional inclusiong of significance values
    
    %make a bar chart and then add error bars without line
    if ~isnumeric(x(1))
        xlabelsticks=x;
        x=1:length(x);
        renamelabels=1;
    else
        renamelabels=0;
    end
    f=figure;
    hold on;
    bar(x,y);
    errorbar(x,y,error_margin,'.','Color','black');

    %then need to add the significance values if they are included
    if nargin>3
        sigs=varargin{1}; %sigs will have combinations of the groups in x e.g. column 1 might be x(1) and column 2 might be x(3), then column 3 has the p value as to whether they are significantly different from one another
        %need to put a star and line between bars with a significance value<0.05
        diffx=abs(sigs(:,1)-sigs(:,2));
        yval=nan(size(sigs,1),1);
        for i =1:size(sigs,1)
        yval(i)=max(y(sigs(i,1):sigs(i,2))+error_margin(sigs(i,1):sigs(i,2)))+3*diffx(i);
        end
        %need to make sure for any y values that are equal that the x values don't overlap, x values are in sigs(,1) and sigs(,2)
        for i=1:size(sigs,1)
            for j=1:size(sigs,1)
                %check if any y values are within +- 3 of each other
                if (abs(yval(i) - yval(j)) <= 3) && i~=j
                    if length((unique([sigs(i,1):sigs(i,2) sigs(j,1):sigs(j,2)])))~=length([sigs(i,1):sigs(i,2) sigs(j,1):sigs(j,2)])
                        if diffx(i)>diffx(j)
                            yval(i)=yval(i)+3;
                        else
                            yval(j)=yval(j)+3;
                        end
                    end
                end
            end
        end

        
        for i=1:size(sigs,1)
            if sigs(i,6)<0.05
                %then need to add a star and line between the two bars
                %get the x values for the two bars
                x1=x(sigs(i,1));
                x2=x(sigs(i,2));
                %get the y values for the two bars
                y1=y(sigs(i,1));
                y2=y(sigs(i,2));
                
                %get the error values for the two bars
                e1=error_margin(sigs(i,1));
                e2=error_margin(sigs(i,2));
                %get the y value for the line
                yline=max(y(sigs(i,1):sigs(i,2))+error_margin(sigs(i,1):sigs(i,2)))+3; %yline=max([y1+e1,y2+e2]);
                %need to make sure the lines are staggered according to the length of bars with longest at the top
                %also need to make sure those of the same length that cross over are staggered oreference the rightmost bar
                lengthb=abs(x1-x2);
                if lengthb>1
                    yline=yline+3*lengthb;
                    
                    
                end
                    yline=yval(i);
                %plot the line and need to make sure it is higher than an
                %of the bars it is crossing
                plot([x1,x2],[yline,yline],'k-');
%                 %plot the veritcal lines leading to the bars
%                 if y1+e1+5>0
%                     plot([x1,x1],[yline,y1+e1+5],'k-');
%                 else
%                     plot([x1,x1],[yline,5],'k-');
%                 end
%                 if y2+e2+5>0
%                     plot([x2,x2],[yline,y2+e2+5],'k-');
%                 else
%                     plot([x2,x2],[yline,5],'k-');
%                 end
                
                %plot the star
                plot(mean([x1,x2]),yline+0.05,'k*');
            end
        end
    end
    %change the x axis labels
    if renamelabels==1
    set(gca,'XTick',1:length(x));
    set(gca,'XTickLabel',xlabelsticks);
    %make sure the x axis labels are at 45 degrees
    set(gca,'XTickLabelRotation',45);
    end
    %set ticks outside
    set(gca,'TickDir','out');
    ylabel('Spike count 4-30ms after stimulus onset');
end