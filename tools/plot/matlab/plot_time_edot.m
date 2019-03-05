function [axis_id] = plot_time_edot (edot, time, varargin)
%--------------------------------------------------------------------------
% [axis_id] = plot_time_edot (edot, time, linecolor)
% Plots the energy generation rate as a function of time, with solid 
% line for positive rate and dotted line for negative rates.  
% Inputs>  edot: time series of energy generation 
%          time: discrete temporal evolution
% Options: LineColor: color of mass fraction lines
%          LineWidth: width of mass fraction lines
%          TimeFromEnd: plot time relative to end of calculation
% Outputs: axis_id: handle of current axis
%--------------------------------------------------------------------------
% Create an instance of the inputParser class.
  p = inputParser;

% Define required inputs
  p.addRequired('edot',@(x)validateattributes(x, {'numeric'}, {'vector', 'real'}));
  p.addRequired('time',@(x)validateattributes(x, {'numeric'}, {'vector', ...
    'real', '>=', 0}));

% Define optional inputs
  p.addOptional(  'LineColor',   'k', @(x)validateattributes(x, {'char'},{'scalar'}));
  p.addOptional(  'LineWidth',     1, @(x)validateattributes(x, {'numeric'}, ...
     {'scalar', 'real', 'positive', '>=', 0.1, '<=', 10}))
  p.addOptional('TimeFromEnd', false, @(x)validateattributes(x, {'logical'},{'scalar'}));
  p.addOptional('LegendOn', true, @(x)validateattributes(x, {'logical'},{'scalar'}));

% Parse and validate all input arguments.
  p.parse(edot, time, varargin{:});

% Assign input arguments
  edot            = p.Results.edot;
  time            = p.Results.time;
  linecolor       = p.Results.LineColor;
  linewidth       = p.Results.LineWidth;
  time_from_end   = p.Results.TimeFromEnd;
  legend_on       = p.Results.LegendOn;

% Build linespecs
  eplus =strcat(linecolor,'-');
  eminus=strcat(linecolor,':');

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
    plot_id=loglog(time,edot,eplus,time,-edot,eminus,'LineWidth',linewidth);

% For large temporal range, us log time axis 
  else      
    plot_id=semilogy(time,edot,eplus,time,-edot,eminus,'LineWidth',linewidth);
    
  end

% Label Plot  
  axis_id = gca;
  ylabel('Energy Production (erg/g/s)');
  xlabel(time_label);
  set(axis_id,'XDir',time_direction);
  set(axis_id,'XLim',[time_start time_stop])
  
% Build Legend
  if(legend_on);
    legend(plot_id,char(' d\epsilon/dt','-d\epsilon/dt'),'Location','NorthEast');
  end
end

