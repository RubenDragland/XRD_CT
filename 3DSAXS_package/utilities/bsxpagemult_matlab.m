rng default

A = double(randi([1 9],3,3,64) );
B = double(randi([1 9],3,8) );
times = 100000;


tic
for ii = times
    R1 = pagemtimes(A,B);
end
toc


tic
for ii = times
    R2 = bsxpagemult(A,B);
end
toc

if any(R1~=R2)
    display(R1);
    display(R2);
end
