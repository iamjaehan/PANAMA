function out = ComputeGini(c)

n = length(c);
sumC = sum(c);

pairSum = 0;
for i = 1:n
    for j = 1:n
        pairSum = pairSum + abs(c(i) - c(j));
    end
end

out = 1/(2*n*sumC)*pairSum;

end