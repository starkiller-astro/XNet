function [] = plot_time_mf_edot( filename, varargin )
%--------------------------------------------------------------------------
%[] = plot_time_mf_edot ( filename, ...)
% Plots mass fractions and energy generation rate vs. time for filename. 
% Inputs>  filename: name of file containing data to plot
% Options: MassLimit: a limiting mass fraction to include in plot
%          MaxSort: when true, species are ordered by maximium mass 
%            fraction; else, lines are ordered as in data file. 
%          TimeFromEnd: plot time relative to end of calculation
%          ZLimit : limits range of species plotted with Z in [Zmin Zmax] 
% Outputs< None
%--------------------------------------------------------------------------
  
% Create an instance of the inputParser class.
  p = inputParser;

% Define required inputs
  p.addRequired('filename',@ischar);

% Define optional inputs
  p.addOptional(  'MassLimit', 1e-25, @(x)validateattributes(x, {'numeric'}, ...
    {'scalar', 'real', 'positive', '>=', 1e-30, '<=', 1e-1}));
  p.addOptional(    'MaxSort',  true, @(x)validateattributes(x, {'logical'},{}));
  p.addOptional('TimeFromEnd', false, @(x)validateattributes(x, {'logical'},{}));
  p.addOptional(     'ZLimit', [0 0], @(x)validateattributes(x, {'numeric'}, ...
                 {'vector', 'integer', 'positive'}));

% Parse and validate all input arguments.
  p.parse(filename, varargin{:});

% Assign input arguments
  filename      = p.Results.filename;
  mass_limit    = p.Results.MassLimit;
  max_sort      = p.Results.MaxSort;
  time_from_end = p.Results.TimeFromEnd;
  z_limit       = p.Results.ZLimit;

% Choose file type
  is_ev_file=strfind(filename,'ev');
  if(isempty(is_ev_file))

% Read TS file
    [zz, aa, xmf, time, ~, ~, ~, edot, ~ , ~] = read_ts_file(filename);

% Build Isotope symbols
    [ nuc_name ] = build_isotope_symbol ( zz,aa );

% Limit abundances plotted by element (zmin < Z <zmax)
    if (z_limit(1)~=0)
      zmin = z_limit(1);
    else
      zmin = min(zz);
    end
    if (z_limit(2)~=0)
      zmax = z_limit(2);
    else
      zmax = max(zz);
    end
    limit_z = find(zz >= zmin & zz<=zmax);
    xmf=xmf(limit_z,:);
    nuc_name=nuc_name(limit_z);

else  
% Alternately, the ev_file may be used
    [ nuc_name, xmf, time, ~, ~, ~, edot] = read_ev_file (filename );

  end

% Limit abundances plotted by maximum mass
  xmf_max= max(xmf,[],2);
  limit_xmf = find(xmf_max > mass_limit);
  xmf=xmf(limit_xmf,:);
  nuc_name=nuc_name(limit_xmf);
  
% Sort isotoped in descending order of maximum mass fraction
  if (max_sort == true)
    xmf_max= max(xmf,[],2);
    [~,sort_order]  = sort(xmf_max,'descend');
    nuc_name=nuc_name(sort_order);
    xmf=xmf(sort_order,:);
  end
  
% Plot time from end if time_from_end is true
  if (time_from_end ==  true)
    time_end = time(size(time,2));
    time=time_end-time;
  end
    
% Plot Mass Fraction
%  line_style='-'
  subplot('position',[.1,.30,.8,.60]);
  mf_axis_id=plot_time_mf ( nuc_name, xmf, time, 'LineWidth', 2, ...
      'TimeFromEnd',time_from_end,'IFile',1);
  set(mf_axis_id,'Ylim',[mass_limit 2])
  set(mf_axis_id,'XTickLabel',[])

% Plot Energy Generation
  subplot('position',[.1,.10,.8,.20]);
  edot_axis_id=plot_time_edot (edot, time, 'LineColor','r','TimeFromEnd',time_from_end);

end

