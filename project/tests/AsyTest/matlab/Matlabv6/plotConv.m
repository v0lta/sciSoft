% Plot the convergence results with LS convergence in loglog, in 1 figure.
% Input
%   r     - Actual values
%   a     - Approximation
%   t     - Title of the plot
%   shift - Degree shift with respect to 2.^(1:maxP2)'
% About
%   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History      - Created October 2013, last edit February 2015
function plotConv(r,a,t,shift)
if ~exist('shift'), shift=0; end
[maxP2,maxOrder] = size(a);
data =  abs(a - repmat(r,1,maxOrder) )./abs(repmat(r,1,maxOrder) );
legends = cell(maxOrder,1);
for i = 1:maxOrder
    if i == 1
        legends{i} = [num2str(i) ' term'];
    else
        legends{i} = [num2str(i) ' terms'];
    end
end
ns = 2.^(1:maxP2)'+shift;
figure; loglog(ns,data,'*'); hold on;
inds = maxP2*ones(1,maxOrder);
tmp = sub2ind(size(data), inds, 1:maxOrder);
% Define `accuracy': which lowest relative error the lines have to interpolate
acc = 10^(-11);
ws = find(data(tmp) < acc );
[rowIdx,colIdx] = find(data(1:maxP2, ws) >= acc );
% with accumarray we take the maximum column index for every row
inds(ws) = accumarray(colIdx,rowIdx,[],@max)';
vals = data(sub2ind(size(data), inds,1:maxOrder) );

loglog(ns, 2.^( (repmat(log2(ns(inds)' ),maxP2,1)-repmat((1:maxP2)',1,maxOrder) ).*...
    (repmat(1:maxOrder,maxP2,1) ) ).*repmat(vals,maxP2,1) );

legend(legends); xlabel('n'); ylabel('Relative error'); title(t);
set(gca, 'XTick', ns);
mir = floor(log10(min(min(data) ) ) );
mar = ceil(log10(max(max(data) ) ) );
set(gca, 'YTick', 10.^(mir:2:mar));
axis([1 2*ns(end) min(min(data))/2 2*max(max(data))]);
