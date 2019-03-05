function [] = compare_time_edot( filenames, varargin )
%--------------------------------------------------------------------------
%[] = compare_time_edot ( filenames, mass_limit, max_sort, ...)
% Plots energy generation rates vs. time for each filename. 
% Inputs>  filenames: a cell array of file names {'file1'; 'file2'; ...}
%          mass_limit: a limiting mass fraction to include in plot
%          max_sort: when true, species are ordered by maximium mass 
%            fraction; else, lines are ordered as in data file. 
% Outputs< None
%--------------------------------------------------------------------------

% Create an instance of the inputParser class.
  p = inputParser;

% Define required inputs
  p.addRequired('filenames',@iscellstr);

% Define optional inputs
  p.addOptional('TimeFromEnd', false, @(x)validateattributes(x, {'logical'},{}));

% Parse and validate all input arguments.
  p.parse(filenames, varargin{:});

% Assign input arguments
  filenames  = p.Results.filenames;
  time_from_end = p.Results.TimeFromEnd
  
% Define manual colors 
  colors=['k','r','g','b','c','m','y'];

  figure;

% Loop over files  
  num_files = size(filenames,1)
%  legend_array=cell(2*num_files)
  
  for ifile= 1:num_files;
    filename = char(filenames(ifile,:));
  
% Choose file type
    is_ev_file=strfind(filename,'ev');
    if(isempty(is_ev_file));
        
% Read TS file
      [~, ~, ~, time, ~, ~, ~, edot, ~ , ~] = read_ts_file(filename);

    else  
% Alternately, the ev_file may be used
      [ ~, ~, time, ~, ~, ~, edot] = read_ev_file (filename );

    end

% Plot Energy Generation
    edot_axis_id=plot_time_edot (edot, time, 'LineColor',colors(ifile),'TimeFromEnd',time_from_end,'LegendOn',false);
    legend_array(2*ifile)  ={strcat(filename,'d\epsilon/dt<0')};
    legend_array(2*ifile-1)={strcat(filename,'d\epsilon/dt>0')};
    hold on
  end
  legend_array;
   [legend_id,hobjs,hlinpatches,legtxt] = legend(edot_axis_id,char(legend_array),'Location','SouthWest');
  hold off
end

