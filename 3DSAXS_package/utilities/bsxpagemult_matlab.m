rng default

A = double(randi([1 9],3,3,3) );
B = double(randi([1 9],3,3,3) );
times = 100000;


tic
for ii = times
    R1 = pagemtimes(A,B);
end
toc
display(R1);

tic
for ii = times
    R2 = bsxpagemult(A,B);
end
toc

display(R2);