function r = isOctave ()
  % https://stackoverflow.com/a/9838357
  persistent x;
  if (isempty (x))
    x = exist ('OCTAVE_VERSION', 'builtin');
  end
  r = x;
end
