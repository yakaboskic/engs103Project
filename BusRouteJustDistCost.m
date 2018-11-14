Cities = ["Hanover" "Lebanon"	"West Lebanon"	"Norwich"	"White River Junction" "Sachem Village"	"Enfield"];
NumStops = 5;
DistCosts = [5.7,4.3,1.6,4.5,1.6,10.8,3.6,7.5,4.5,7,6.7,5.6,1,3.6,10.7,5.2,3.2,12.8,3.8,11.5,12.3];
lendist = length(DistCosts);
Pops = [0, 98, 42, 26,	36,	72,	13];
f = DistCosts;
idxs = nchoosek(1:7,2);

Aeq = spones(1:length(idxs));
beq = NumStops;

Aeq = [Aeq;spalloc(NumStops,length(idxs),NumStops*(NumStops-1))]; 
for ii = 1:NumStops
    whichIdxs = (idxs == ii); 
    whichIdxs = sparse(sum(whichIdxs,2)); 
    Aeq(ii+1,:) = whichIdxs'; 
end
beq = [beq; 2*ones(NumStops,1)];

lb = zeros(length(f),1);
ub = ones(length(f),1);
intcon = 1:length(f);

[x,fval,exitflag,output] = intlinprog(f,intcon,[],[],Aeq,beq,lb,ub);

tours = detectSubtours(x,idxs);
numtours = length(tours);

Asub = spalloc(0,length(f),0); 
bsub = [];
while numtours > 1 
    bsub = [bsub;zeros(numtours,1)]; 
    A1 = spalloc(numtours,lendist,NumStops);
    Asub = [Asub;A1]; 
    for ii = 1:numtours
        rowIdx = size(Asub,1)+1;
        subTourIdx = tours{ii}; 
        variations = nchoosek(1:length(subTourIdx),2);
        for jj = 1:length(variations)
            whichVar = (sum(idxs==subTourIdx(variations(jj,1)),2)) & ...
                       (sum(idxs==subTourIdx(variations(jj,2)),2));
            Asub(rowIdx,whichVar) = 1;
        end
        bsub(rowIdx) = length(subTourIdx)-1;
    end
    [x,fval,exitflag,output] = intlinprog(f,intcon,Asub,bsub,Aeq,beq,lb,ub);

    tours = detectSubtours(x,idxs);
    numtours = length(tours); 
end
Cities(idxs(x==1,:))
Stops = Cities(unique(idxs(x==1,:)))
PopMet = sum(Pops(unique(idxs(x==1,:))))
fval