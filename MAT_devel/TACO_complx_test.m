nSet = [2, 5, 10, 20, 30, 40, 50];
mSet = [10, 50, 100, 500];
gamma = 0.9;
epsilon = 1e-3;

timeRecord = zeros(length(nSet), length(mSet));
iteration = 5;

for i = nSet
    for j = mSet
        for ii = 1:iteration
        n = i;
        m = j;
        C = rand(n,m);
        b = rand(n,1);
        d0 = 1;
            disp([num2str(i)+", "+num2str(j)])
            tic
            TACo(C, b, d0, gamma, epsilon);
            time = toc;
            timeRecord(find(nSet==i),find(mSet==j)) = timeRecord(find(nSet==i),find(mSet==j)) + time;
        end
    end
end

timeRecord = timeRecord./iteration;

figure(1); clf;
mesh(nSet, mSet,log10(timeRecord))
xlabel("n")
ylabel("m")
zlabel("logtime [s]")
