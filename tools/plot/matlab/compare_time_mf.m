function [] = compare_time_mf( filenames, varargin )
%--------------------------------------------------------------------------
%[] = compare_time_mf ( filenames, mass_limit, max_sort, ...)
% Plots mass fractions vs. time for each filename. 
% Inputs>  filenames: a cell array of file names {'file1'; 'file2'; ...}
%          MassLimit: a limiting mass fraction to include in plot
%          MaxSort: when true, species are ordered by maximium mass 
%            fraction; else, lines are ordered as in data file. 
%          ZLimit: limit range in proton number to be graphed
%          TimeFromEnd: plot time relative to end of calculation
%          SeparateLegend: draw separate legends for each lineset
% Outputs< None
%--------------------------------------------------------------------------

% Create an instance of the inputParser class.
  p = inputParser;

% Define required inputs
  p.addRequired('filenames',@iscellstr);

% Define optional inputs
  p.addOptional(  'MassLimit', 1e-25, @(x)validateattributes(x, {'numeric'}, ...
                {'scalar', 'real', 'positive', '>=', 1e-30, '<=', 1e-1}));
  p.addOptional(       'MaxSort',  false, @(x)validateattributes(x, {'logical'},{}));
  p.addOptional(   'TimeFromEnd', false, @(x)validateattributes(x, {'logical'},{}));
  p.addOptional('SeparateLegend', false, @(x)validateattributes(x, {'logical'},{}));
  p.addOptional(     'ZLimit', [0 0], @(x)validateattributes(x, {'numeric'}, ...
    {'vector', 'integer'}));

% Parse and validate all input arguments.
  p.parse(filenames, varargin{:});

% Assign input arguments
  filenames       = p.Results.filenames;
  mass_limit      = p.Results.MassLimit;
  max_sort        = p.Results.MaxSort;
  time_from_end   = p.Results.TimeFromEnd;
  separate_legend = p.Results.SeparateLegend;
  z_limit         = p.Results.ZLimit;
  
% Loop over files  
  num_files = size(filenames,1);
  for ifile= 1:num_files;
    filename = char(filenames(ifile,:));
  
% Choose file type
    is_ev_file=strfind(filename,'ev');
    if(isempty(is_ev_file));
        
% Read TS file
      [zz, aa, xmf, time, ~, ~, ~, ~, ~ , ~] = read_ts_file(filename);

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
      [ nuc_name, xmf, time, ~, ~, ~, ~] = read_ev_file (filename );

    end

% Limit abundances plotted to ncolor largest in maximum mass
    ncolor = 14;
    xmf_max= max(xmf,[],2);
    [~,sort_order]  = sort(xmf_max,'descend');
% Sort isotopes in descending order of maximum mass fraction
    if (max_sort == true);
      max_value_indices = sort_order(1:ncolor);
    else
      max_value_indices = sort(sort_order(1:ncolor),'ascend');
    end
    xmf=xmf(max_value_indices,:);
    nuc_name=nuc_name(max_value_indices);
  
% Warn if more species than we can compare
    desired_species = size(sort_order,1);
    if (desired_species > ncolor);
        warning('Number of desired species, %d, too large, truncating to number of colors, %d!', desired_species, ncolor);
    end
        
% Plot Mass Fraction
    mf_axis_id=plot_time_mf (nuc_name, xmf, time, 'LineWidth', 2, ...
      'TimeFromEnd',time_from_end,'IFile',ifile,'SeparateLegend',separate_legend); 
    set(mf_axis_id,'Ylim',[mass_limit 2])
    hold on
  end
  hold off
end

