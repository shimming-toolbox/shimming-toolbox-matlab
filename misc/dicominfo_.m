function hdr = dicominfo_(varargin)
  % for some reason, the octave version of dicominfo is giving a single trailing ' ' on many (all?) fields.
  hdr = dicominfo(varargin{:});
  fields = fieldnames(hdr);
  for i=1:numel(fields)
    if isa(hdr.(fields{i}), 'char')
      hdr.(fields{i}) = strtrim(hdr.(fields{i}));
    end
  end

end
