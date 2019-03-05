function [ nuc_name ] = build_isotope_symbol ( zz,aa )
%--------------------------------------------------------------------------
% [nuc_name] = build_isotope_symbol ( zz,aa ) 
% Returns TeX string symbols for each isotope including superscripting.
% Inputs>  zz: zz: proton number for each isotope in the set
%          aa: total nucleon number for each isotope in the set
% Outputs< nuc_name: cell array of isotopic symbols
%--------------------------------------------------------------------------

% Get element symbols
  [element] = build_element_symbol;

% Build the isotope symbols
  n_iso=size(zz,1);
  for i=1:n_iso
    [label,errmsg]=sprintf('^{%d}%s',aa(i),element(zz(i)+1,:));
     nuc_name(i)={label};
  end

end

