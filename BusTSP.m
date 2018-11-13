Cities = ["Hanover" "Lebanon"	"West Lebanon"	"Norwich"	"White River Junction" "Sachem Village"	"Enfield"];
idxs = nchoosek(1:7,2);
nStops = 7;
D = [30	5.7	4.3	1.6	4.5	1.6	10.8
5.7	30	3.6	7.5	4.5	7	6.7
4.3	3.6	30	5.6	1	3.6	10.7
1.6	7.5	5.6	30	5.2	3.2	12.8
4.5	4.5	1	5.2	30	3.8	11.5
1.6	7	3.6	3.2	3.8	30	12.3
10.8	6.7	10.7	12.8	11.5	12.3	30];
dist = [];
for i = 1:length(idxs)
    dist = [dist,D(idxs(i,1),idxs(i,2))];
end
lendist = length(dist);

Aeq = spones(1:lendist); % Adds up the number of trips
beq = nStops;

Aeq = [Aeq;spalloc(nStops,length(idxs),nStops*(nStops-1))]; % allocate a sparse matrix
for ii = 1:nStops
    whichIdxs = (idxs == ii); % find the trips that include stop ii
    whichIdxs = sparse(sum(whichIdxs,2)); % include trips where ii is at either end
    Aeq(ii+1,:) = whichIdxs'; % include in the constraint matrix
end
beq = [beq; 2*ones(nStops,1)];
intcon = 1:lendist;
lb = zeros(lendist,1);
ub = ones(lendist,1);

opts = optimoptions('intlinprog','Display','off','Heuristics','round-diving',...
    'IPPreprocess','none');
[x_tsp,costopt,exitflag,output] = intlinprog(dist,intcon,[],[],Aeq,beq,lb,ub,opts);

Cities(idxs(x_tsp==1,:))
