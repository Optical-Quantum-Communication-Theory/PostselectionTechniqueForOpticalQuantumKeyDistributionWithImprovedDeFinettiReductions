load temp.mat



eqConstViolations = zeros(numel(eqCons),1);

for index=1:numel(eqCons)
    oper = eqCons(index).operator;
    val = eqCons(index).scalar;
    eqConstViolations(index) = abs(trace(oper*rho) - val);
end

ineqConstViolations = zeros(numel(ineqCons),1);

for index=1:numel(eqCons)
    oper = ineqCons(index).operator;
    lowerval = ineqCons(index).lowerBound;
    upperval = ineqCons(index).upperBound;

    ineqConstViolations(index) = max( max((trace(oper*rho) - upperval),0),...
        max(lowerval - trace(oper*rho),0));
end
