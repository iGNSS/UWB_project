function [flag,fifo] = buf(data,fifo,vol)
%FIFO；
[row,~] = size(fifo);
flag = 0;
if row < vol
    fifo(row+1,:) = data;
else
    fifo(1,:) = [];
    fifo(row,:) = data;
end
if row >= (vol-1)
    flag = 1;%FIFO存满；
end
end