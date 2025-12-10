nSet = [3, 10, 30, 100];
mSet = [3, 10, 30, 100, 300, 1000];
% nSet = 3;
% mSet = 3;
gamma = 0.9;
epsilon = 1e-2;
iteration = 100;
d0 = 1;

timeRecord = zeros(length(nSet), length(mSet), iteration);
roundRecord = zeros(length(nSet), length(mSet), iteration);
cycleStep = zeros(length(nSet), length(mSet), iteration);

for i = nSet
    for j = mSet
        for ii = 1:iteration
            n = i;
            m = j;
            C = rand(n,m);
            b = rand(n,1);
            if mod(ii,10) == 0
                disp([num2str(i)+", "+num2str(j)+", "+num2str(ii)])
            end
            tic
            [~,~,~,rounds,maxCycleStep] = TACo(C, b, d0, gamma, epsilon);
            time = toc;
            % timeRecord(find(nSet==i),find(mSet==j),ii) = timeRecord(find(nSet==i),find(mSet==j)) + time;
            timeRecord(find(nSet==i),find(mSet==j),ii) = time;
            roundRecord(find(nSet==i),find(mSet==j),ii) = rounds;
            cycleStep(find(nSet==i),find(mSet==j),ii) = maxCycleStep;
        end
    end
end

% timeRecord = timeRecord./iteration;

% figure(1); clf;
% mesh(nSet, mSet,log10(timeRecord))
% xlabel("n")
% ylabel("m")
% zlabel("logtime [s]")
%%
markers = {'o','s','^','d'};  % 각 라인별 마커 지정
figure(2)
clf
for i = 1
    localInfo = roundRecord(i,:,:);
    % localInfo = timeRecord(i,:,:);
    localCycleInfo = cycleStep(i,:,:);
    errorbar(mSet, mean(localInfo,3), std(localInfo,[],3)/sqrt(iteration),'LineWidth',2,'Marker','.','MarkerSize',20)
    hold on
    % plot(mSet, max(localInfo,[],3),'LineWidth',2,'Marker','+','LineStyle',':')
    % plot(mSet, max(localCycleInfo,[],3),'LineWidth',2,'Marker','+','LineStyle',':')
end
grid on
xlabel("Number of choices ($m$)",'Interpreter','latex')
ylabel("Number of rounds",'Interpreter','latex')
ylim([1 inf])
set(gca,'YScale','log','XScale','log','FontName','times','fontsize',23)
set(gcf,'Position',[100 100 800 600])
legend("$n = 3$","$n = 10$","$n = 30$","$n = 100$","Location","northwest",'interpreter','latex');
set(gca,'FontName','times','FontSize',23)
xlim([2 1300])
% ylim([1 8000])
exportgraphics(gca,'./tacoComplex.pdf','Resolution',300)

% figure(3)
% clf
% for i = 1:3
%     localInfo = roundRecord(i,:,:);
%     % localInfo = timeRecord(i,:,:);
%     errorbar(mSet, mean(localInfo,3), std(localInfo,[],3),'LineWidth',3)
%     hold on
% end
% grid on
% ub = zeros(1,length(mSet));
% ub_temp = zeros(1,length(mSet));
% for i = 1:length(mSet)
%     m = mSet(i);
%     n = 3;
%     C = 3;
%     ub_temp(i) = log(1/(m+1)/d0/(n-1))/log(gamma);
%     % ub(i) = ub_temp(i) * (m+1) * n * m / 10;
%     % ub(i) = ub_temp(i) * (n+1) * n * n / 10000;
%     % ub(i) = (m+1) * n * m;
%     % ub(i) = ub_temp(i)*1;
%     ub(i) = ub_temp(i) * 1;
%     for j = 1:n-1
%         m = n;
%         ub(i) = ub(i) * nchoosek(C+m-n+j-1, m-1);
%     end
% end
% plot(mSet,ub,'r--')
% set(gca,'YScale','log','XScale','log')