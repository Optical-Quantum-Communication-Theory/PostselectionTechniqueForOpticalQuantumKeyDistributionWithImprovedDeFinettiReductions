function PSPlots(liftType)


Colour = {"#0072BD"	,"#D95319" ,"#EDB120", "#7E2F8E","#77AC30",	"#4DBEEE", "#A2142F"};
styles = {'-', '--', ':', '-.','-', '--', ':'}; % Line styles
markers = {'o', 'x', '*', '+', '.','square', 'diamond'}; % Marker types

for i=1:numel(liftType)
    data = load("FiniteThreeState"+liftType{i}+".mat");
    keyRate = [data.results.keyRate];
    keyRate(keyRate<0) = 0;
    for iter = 1:numel(keyRate)
        keyLength(iter) = [data.results(iter).debugInfo.keyRateModule.keyLength];
    end
    keyLength(keyLength<0) = 0;

    distanceList = -10*log10(cell2mat(data.qkdInput.scanParameters.eta))/0.16;

    semilogy(distanceList ,keyLength, [markers{i},styles{i}], "Color",Colour{i}, 'DisplayName',liftType{i});
    hold on;

end


%Other things

fsize = 42;
lw = 4;
ms = 13;


xlabel('distance (km)','Interpreter','latex');
ylabel('key length','Interpreter','latex');
legend('Location','best','Interpreter','latex');
fontsize(gca,fsize,"pixels")

set(findall(gcf,'Type','line'),'LineWidth',lw);
set(findall(gcf,'Type','line'),'MarkerSize',ms);

hold off;
end
