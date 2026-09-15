% ----------------------------------------------------------------------- %
% Function plot_areaerrorbar_multiple plots the mean and standard         %
% deviation of a set of data filling the space between the positive and   %
% negative mean error using a semi-transparent background, completely     %
% customizable.                                                           %
%                                                                         %
%   Input parameters:                                                     %
%       - **data:   Cell of data matrices, with rows corresponding to     %
%                   observations and columns to samples.                  %
%       - options:  (Optional) Struct that contains the customized params.%
%           * options.handle:       Figure handle to plot the result.     %
%           ** options.legend:      Name of legend entries.               %
%           * options.color_area:   RGB color of the filled area.         %
%           * options.color_line:   RGB color of the mean line.           %
%           * options.alpha:        Alpha value for transparency.         %
%           * options.line_width:   Mean line width.                      %
%           * options.x_axis:       X time vector.                        %
%           * options.error:        Type of error to plot (+/-).          %
%                   if 'std',       one standard deviation;               %
%                   if 'sem',       standard error mean;                  %
%                   if 'var',       one variance;                         %
%                   if 'c95',       95% confidence interval.              %
% ----------------------------------------------------------------------- %
%   Example of use:                                                       %
%       data = repmat(sin(1:0.01:2*pi),100,1);                            %
%       data = data + randn(size(data));                                  %
%       plot_areaerrorbar(data);                                          %
% ----------------------------------------------------------------------- %
%   Author:  Victor Martinez-Cagigal                                      %
%   Date:    30/04/2018                                                   %
%   E-mail:  vicmarcag (at) gmail (dot) com                               %
%   Modified by Mia Brunetti                                              %
%   Date:    02/10/2021                                                   %
%   Modification: can add multiple graphs on the same plot and displays a %
%   legend of only the means.                                             %
% ----------------------------------------------------------------------- %
function plot_areaerrorbar_multiple(data_cell, options)

[row,col] = cellfun(@size, data_cell);
data_mean = zeros(size(data_cell,1),col(1));
data_std = zeros(size(data_cell,1),col(1));
data_error = zeros(size(data_cell,1),col(1));

    for i = 1:size(data_cell,1)
        data = data_cell{i,1};

        if(isfield(options,'x_axis')==0)
            options.x_axis = 1:size(data,2); 
        end
        options.x_axis = options.x_axis(:);

        % Computing the mean and standard deviation of the data matrix
        data_mean(i,:) = mean(data,1);
        data_std(i,:)  = std(data,0,1);

        % Type of error plot
        switch(options.error)
            case 'std', data_error(i,:) = data_std(i,:);
            case 'sem', data_error(i,:) = (data_std(i,:)./sqrt(size(data,1)));
            case 'var', data_error(i,:) = (data_std(i,:).^2);
            case 'c95', data_error(i,:) = (data_std(i,:)./sqrt(size(data,1))).*1.96;
        end

        % Plotting the result
        x_vector(i,:) = [options.x_axis(1:10:end)', fliplr(options.x_axis(1:10:end)')];
        y_vector(i,:) = [data_mean(i,1:10:end)+data_error(i,1:10:end),fliplr(data_mean(i,1:10:end)-data_error(i,1:10:end))];
    
    end
    
    % Plotting the result
    color_area = options.color_area;
    color_line = options.color_line;
    doses = options.legend;
    
    %figure(options.handle);
    hold on;
    for j = 1:size(data_cell,1)
        shape = fill(x_vector(j,:),y_vector(j,:), color_area(j,:), 'DisplayName', '');
        shape.Annotation.LegendInformation.IconDisplayStyle = 'off';
        set(shape, 'edgecolor', 'none');
        set(shape, 'FaceAlpha', options.alpha);
        txt = [sprintf('%g nM',doses(1,j))]; %options.legend inputs are initial values (doubles)
        plot(options.x_axis, data_mean(j,:), 'color', color_line(j,:), ...
            'LineWidth', options.line_width, 'DisplayName', txt);
    end
    hold off;
end