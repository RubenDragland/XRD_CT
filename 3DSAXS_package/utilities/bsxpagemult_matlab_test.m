

a = ones(3,3,3);
Ylm = ones(3,3,3);

tic
for ii = 100000
    out = bsxpagemult(a, Ylm);
end
toc

fmt=['out =' repmat(' %1.0f',1,numel(out))];
fprintf(fmt, out);

