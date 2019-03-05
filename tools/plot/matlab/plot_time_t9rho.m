function [axis_id] = plot_time_t9rho (temperature, density, time, varargin)
%--------------------------------------------------------------------------
% [axis_id] = plot_time_t9rho (temperature, density, time, linecolor)
% Plots the energy generation rate as a function of time, with solid ...
% line for positive rate and dotted line for negative rates.  
% Inputs>  temperature: time series of temperature 
%          density: time series of density 
%          time: discrete temporal evolution
% Options: LineColor: color of mass fraction lines
%          LineWidth: width of mass fraction lines
%          TimeFromEnd: plot time relative to end of calculation
% Outputs< axis_id: handle of current axis
%--------------------------------------------------------------------------

% Create an instance of the inputParser class.
  p = inputParser;

% Define required inputs
  p.addRequired('temperature', @(x)validateattributes(x, {'numeric'}, {'vector', 'real'})); 
  p.addRequired(    'density', @(x)validateattributes(x, {'numeric'}, {'vector', 'real'}));
  p.addRequired(       'time', @(x)validateattributes(x, {'numeric'}, {'vector', 'real', '>=', 0}));

% Define optional inputs
  p.addOptional(  'LineColor',   'k', @(x)validateattributes(x, {'char'},{'scalar'}));
  p.addOptional(  'LineWidth',     1, @(x)validateattributes(x, {'numeric'}, ...
               {'scalar', 'real', 'positive', '>=', 0.1, '<=', 10}))
  p.addOptional('TimeFromEnd', false, @(x)validateattributes(x, {'logical'},{'scalar'}));

% Parse and validate all input arguments.
  p.parse(temperature, density, time, varargin{:});

% Assign input arguments
  temperature     = p.Results.temperature;
  density         = p.Results.density;
  time            = p.Results.time;
  linecolor       = p.Results.LineColor;
  linewidth       = p.Results.LineWidth;
  time_from_end   = p.Results.TimeFromEnd;

% Test dynamic range in time
  ntime=size(time,2);
  if( time_from_end == true);
    time_stop=time(1);
    time_start=time(ntime-1);
    time_range=time_start/time_stop;
    time_label='Time from Completion (s)';
    time_direction='reverse';
  else
    time_start=time(2);
    time_stop=time(ntime);
    time_range=time_stop/time_start;
    time_label='Time(s)';
    time_direction='normal';
  end

% For small temporal range, us linear time axis 
  if(time_range > 20) 
    [axis_id,line1_id,line2_id] = plotyy(time,temperature,time,density,'semilogx','loglog');

% For large temporal range, us log time axis 
  else      
    [axis_id,line1_id,line2_id]=plotyy(time,temperature,time,density,'plot','semilogy');
    
  end
  
  % Label Plot
  xlabel(time_label);
  set(axis_id(1),'XDir',time_direction);
  set(axis_id(1),'XLim',[time_start time_stop]);
  set(axis_id(1),'YTickMode','auto');
  set(axis_id(1),'YMinorTick','on');
  set(axis_id(2),'XDir',time_direction);
  set(axis_id(2),'XLim',[time_start time_stop]);
  set(axis_id(2),'YTickMode','auto');
  set(get(axis_id(1),'Ylabel'),'String','Temperature (GK)');
  set(get(axis_id(2),'Ylabel'),'String','Density (g/cm^3)');
  set(axis_id(2),'YColor','k');

% Set linespecs
  set(line1_id,'LineStyle','-');
  set(line1_id,'LineWidth',linewidth);
  set(line1_id,'Color',linecolor);
  set(line2_id,'LineStyle','--');
  set(line2_id,'LineWidth',linewidth);
  set(line2_id,'Color',linecolor);

% Build Legend
  legend(axis_id(1),char('Temperature','Density    '),'Location','NorthEast');

end

