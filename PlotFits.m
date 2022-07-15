function PlotFits(data, fits, params)
% Plots the fits produced by the dislocation fitting package, according to
% the params specified
%

for i=1:length(fits) %iterate through all the fits
    
    
    if isempty(fits(i).gof) %only plot ones with goodness of fit data
        continue; 
    end

    tempgof = fits(i).gof; %retrieve gof data for fit i; kind of a janky work around to access data 
                                %stored in a substructure
            
    figure;
    hold on;

    %Plot all the data with error bars
    if params.ErrBars
        ebs1 = errorbar(data.angles, data.values, error1, '.b');
%             errorbar(load2.data.angles, load2.data.values, error2, '.b');
%             errorbar(load3.data.angles, load3.data.values, error3, '.b');

    else
        ebs1 = plot(data.angles, data.values, '.b');
%             plot(load2.data.angles, load2.data.values, '.b');
%             plot(load3.data.angles, load3.data.values, '.b');
    end

    %Plot the fit
    fithandle = plot(fits(i).fits);
    fithandle.LineWidth = 1.5;

    xlabel('\theta (rad)','FontSize',18,'FontWeight','bold');
    ylabel('Displacement (ang.)','FontSize',18,'FontWeight','bold');
    prop = gca;
    prop.FontWeight = 'bold';
    prop.FontSize = 18;
    legend([ebs1,fithandle], 'Data', 'Fit', 'Location', 'southeast');
    legend('boxoff');
    ylim( params.YRange )

    %Get length of Burger's vector for reference
    blength = norm(ObjectiveCoordTransform(params.lat, fits(i).ijk))/params.prefactor;

    %If short boi save figure to a file before annotating
    if(blength - params.Burgers_ciel)<0.1
        filename = [params.UnprocDir '\fit2 ' num2str(fits(i).ijk) '.png'];
        saveas(gcf, filename); 
    end

    %Add some formatting for the figure
    title(['[' num2str(fits(i).ijk) ']']);

    %Annotate plot with coefficients and confidence intervals
    if params.Annotate
        coeffs = coeffnames(fits(i).fits);
        coeffvals= coeffvalues(fits(i).fits);
        ci = confint(fits(i).fits,0.95);
        str1 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{1},coeffvals(1),ci(:,1));
        str2 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{2},coeffvals(2),ci(:,2));
        str3 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{3},coeffvals(3),ci(:,3));
        str4 = sprintf('\n Burgers Vector Length = %0.3f', blength);
        str5 = sprintf('\n RMSE of Fit = %0.3f', tempgof.rmse);
        annotation('textbox',[.15 .5 .7 .4],'String',...
            ['Coefficients (with 95% confidence bounds): ', str1, str2, str3, str4, str5], 'EdgeColor', 'none');
    %         fits(i).ijk
    %         fits(i).gof
    end
    hold off;

    %Close all the ones with long Burgers vectors
    if (params.Burgers_ciel-blength)<0; close; end

    continue       
end


end

