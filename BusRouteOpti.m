%Cities = ["Hanover" "Lebanon"	"West Lebanon"	"Norwich"	"White River Junction" "Sachem Village"	"Enfield"];
alpha = 1;
beta = 1;
NumStops = 7;
DistCosts = [5.7,4.3,1.6,4.5,1.6,10.8,3.6,7.5,4.5,7,6.7,5.6,1,3.6,10.7,5.2,3.2,12.8,3.8,11.5,12.3];
lendist = length(DistCosts);
Pops = [0, 98,	42,	26,	36,	72,	13];
f = cat(2,alpha*DistCosts,-beta*Pops);
idxs = nchoosek(1:7,2);

A = [0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	1	1	1	1	1
-1	-1	-1	-1	-1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0
-1	0	0	0	0	0	-1	-1	-1	-1	-1	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0
0	-1	0	0	0	0	-1	0	0	0	0	-1	-1	-1	-1	0	0	0	0	0	0	0	0	1	0	0	0	0
0	0	-1	0	0	0	0	-1	0	0	0	-1	0	0	0	-1	-1	-1	0	0	0	0	0	0	1	0	0	0
0	0	0	-1	0	0	0	0	-1	0	0	0	-1	0	0	-1	0	0	-1	-1	0	0	0	0	0	1	0	0
0	0	0	0	-1	0	0	0	0	-1	0	0	0	-1	0	0	-1	0	-1	0	-1	0	0	0	0	0	1	0
0	0	0	0	0	-1	0	0	0	0	-1	0	0	0	-1	0	0	-1	0	-1	-1	0	0	0	0	0	0	1];

b = zeros(size(A,1),1);
b(1) = NumStops;

Aeq = [0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0
1	1	1	1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
1	0	0	0	0	0	1	1	1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	1	0	0	0	0	1	0	0	0	0	1	1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	1	0	0	0	0	1	0	0	0	1	0	0	0	1	1	1	0	0	0	0	0	0	0	0	0	0
0	0	0	1	0	0	0	0	1	0	0	0	1	0	0	1	0	0	1	1	0	0	0	0	0	0	0	0
0	0	0	0	1	0	0	0	0	1	0	0	0	1	0	0	1	0	1	0	1	0	0	0	0	0	0	0
0	0	0	0	0	1	0	0	0	0	1	0	0	0	1	0	0	1	0	1	1	0	0	0	0	0	0	0];

beq = 2*ones(size(A,1),1);
beq(1) = 1;
lb = zeros(28,1);
ub = ones(28,1);
intcon = [1:numel(f)];

[x,fval,exitflag,output] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub);

tours = detectSubtours(x(1:21),idxs);
numtours = length(tours);

A = spalloc(0,lendist,0); % Allocate a sparse linear inequality constraint matrix
b = [];
while numtours > 1 % repeat until there is just one subtour
    % Add the subtour constraints
    b = [b;zeros(numtours,1)]; % allocate b
    A = [A;spalloc(numtours,length(f),NumStops)]; % a guess at how many nonzeros to allocate
    for ii = 1:numtours
        rowIdx = size(A,1)+1; % Counter for indexing
        subTourIdx = tours{ii}; % Extract the current subtour
%         The next lines find all of the variables associated with the
%         particular subtour, then add an inequality constraint to prohibit
%         that subtour and all subtours that use those stops.
        variations = nchoosek(1:length(subTourIdx),2);
        for jj = 1:length(variations)
            whichVar = (sum(idxs==subTourIdx(variations(jj,1)),2)) & ...
                       (sum(idxs==subTourIdx(variations(jj,2)),2));
            A(rowIdx,whichVar) = 1;
        end
        b(rowIdx) = length(subTourIdx)-1; % One less trip than subtour stops
    end

    % Try to optimize again
    [x,costopt,exitflag,output] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub);
    
    % How many subtours this time?
    tours = detectSubtours(x(1:21),idxs);
    numtours = length(tours); % number of subtours
    fprintf('# of subtours: %d\n',numtours);
end
x
costopt
