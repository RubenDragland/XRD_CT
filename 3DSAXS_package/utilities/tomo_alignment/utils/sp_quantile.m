function Q = sp_quantile(x,q,reduce)
   % sparse quantile (make a fast guess for large arrays)
   x = x(1:reduce:end);
   Q = quantile(x,q); 
end