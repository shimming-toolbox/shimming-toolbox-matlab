function [val, new_fp] = mread(arr, fp, count, type)
% mread: eja (almost) drop-in replacement for fread for reading uint8 arrays
% dependencies: c_str.m

[computerType, maxSize, endian] = computer; %#ok<ASGLU>
isLittleEndian = (endian == 'L');

if (strcmpi(type, 'int32-le'))
    % data type is int32, little endian
    val = zeros(count, 1, 'int32');
    for x=1:count
        cur_ptr = fp + (x-1)*4;
        tmp_dat = arr(cur_ptr:cur_ptr+4-1);
        if (isLittleEndian)
            val(x) = typecast(tmp_dat, 'int32');
        else
            val(x) = swapbytes(typecast(tmp_dat, 'int32'));
        end
    end
    val = double(val);
    new_fp = fp + count*4;
elseif (strcmpi(type, 'uint32-le'))
    % data type is uint32, little endian
    val = zeros(count, 1, 'uint32');
    for x=1:count
        cur_ptr = fp + (x-1)*4;
        tmp_dat = arr(cur_ptr:cur_ptr+4-1);
        if (isLittleEndian)
            val(x) = typecast(tmp_dat, 'uint32');
        else
            val(x) = swapbytes(typecast(tmp_dat, 'uint32'));
        end
    end
    val = double(val);
    new_fp = fp + count*4;
elseif (strcmpi(type, 'c_str'))
    % data type is null-terminated string
    val = char(arr(fp:fp+count-1));
    if (size(val,1) > size(val,2)), val = val'; end
    val = c_str(val);
    new_fp = fp + count;
elseif (strcmpi(type, 'char'))
    % data type is char
    val = char(arr(fp:fp+count-1));
    if (size(val,1) > size(val,2)), val = val'; end
    new_fp = fp + count;
else
    error('mread: unknown data type!');
end
