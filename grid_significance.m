function grid_significance(x,significanceval)
%creates an X x X grid of significance values with stars in each square to indicate the level of significance
%non-significant values are left blank
%significant values are also shaded according to the level of significance
%significanceval has the x values in the first and second column that it is comparing e.g. condition 1 of x against condition 7 would put 1 in the first column and 7 in the second column
%the sixth column has the significance value
grid=nan(length(x),length(x));
for i=1:length(x)
    for j=1:length(x)
        if i~=j
            for k=1:size(significanceval,1)
                if (significanceval(k,1)==i && significanceval(k,2)==j) ||  (significanceval(k,1)==j && significanceval(k,2)==i)
                    grid(i,j)=significanceval(k,6);
                end
            end
        end
    end
end



%plot the grid
figure;
h=imagesc(grid); % Apply logarithmic transformation to grid values
colormap(parula); % Use a linear colormap
colorbar;
set(h, 'AlphaData', ~isnan(grid))




%set the x and y labels
set(gca,'XTick',1:length(x));
set(gca,'YTick',1:length(x));
set(gca,'XTickLabel',x);
set(gca,'YTickLabel',x);
%set the x and y labels at 45 degrees
set(gca,'XTickLabelRotation',45);
hold on;
%plot the stars in each square of the grid according to significance 
for i=1:length(x)
    for j=1:length(x)
        if i~=j
            if grid(i,j)<0.001
                text(j,i,'***','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',15);
            elseif grid(i,j)<0.01
                text(j,i,'**','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',15);
            elseif grid(i,j)<0.05
                text(j,i,'*','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',15);
            end
        end
    end
end
caxis([0 0.05]); % Adjust the color axis limits

set(gca,'ColorScale','log')
%set the axis to be equal
axis equal;
%set the axis to be tight
axis tight;

%set the ticks to be out
set(gca,'TickDir','out');

end