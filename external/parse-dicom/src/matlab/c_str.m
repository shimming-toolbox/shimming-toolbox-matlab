function outstr = c_str(instr)
% eja function to read C-style null-terminated strings
% (strips everything past the first null)

outstr = strtok(instr, char(0));
