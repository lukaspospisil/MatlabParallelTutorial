clear all
close all

% prepare grid
h = 2e-2;
[X,Y] = meshgrid(0:h:10,0:h:10);

disp(['Number of evaluations (work to do): ' num2str(size(X,1)*size(X,2))])

%% SEQUENTIAL
FXY_Sequential = zeros(size(X));
tStartSequential = tic;
disp('-> solving sequential problem')
for i=1:size(FXY_Sequential,1)
    for j=1:size(FXY_Sequential,2)
        FXY_Sequential(i,j) = evaluate_something(X(i,j),Y(i,j));
    end
end
timeSequential = toc(tStartSequential);

%% PARALLEL
numOfWorkers = 8; % number of parallel processes
poolobj = parpool(numOfWorkers); % start parallel environment

% unfortunatelly, it is not possible to nest parfors,
% therefore we have to shrink them into one parfor
% idea: separate data into disjoint parts
% every process will operate with one {i}

% matrix to vector (for easier separation)
X_size = size(X,1); % number of rows
Y_size = size(X,2); % number of columns
globalsize = X_size*Y_size; % the amount of all work which has to be computed
Xvec = reshape(X,[globalsize,1]); % matrix X into vector
Yvec = reshape(Y,[globalsize,1]); % matrix Y into vector

% local variables for each worker 
idxparallel = cell(numOfWorkers,1); % indexes of the work for every worker
Xparallel = cell(numOfWorkers,1); % worker portion of X
Yparallel = cell(numOfWorkers,1); % worker portion of Y
FXYparallel = cell(numOfWorkers,1); % worker portion of FXY_Parallel

% estimate size of local portion of data owned by each worker (optimal)
localsize = ceil(globalsize/numOfWorkers);

% separate data
for proc=1:numOfWorkers
    idxparallel{proc} = (proc-1)*localsize+1:min(proc*localsize,globalsize); % compute local indexes
    Xparallel{proc} = Xvec(idxparallel{proc});
    Yparallel{proc} = Yvec(idxparallel{proc});
    FXYparallel{proc} = zeros(length(idxparallel{proc}),1); % initialize
end

% run in parallel
tStartParallel = tic;
disp('-> solving parallel problem')
parfor proc=1:numOfWorkers
    % process local portion of work
    for idxlocal=1:length(idxparallel{proc});
        FXYparallel{proc}(idxlocal) = evaluate_something(Xparallel{proc}(idxlocal),Yparallel{proc}(idxlocal));
    end
end
timeParallel = toc(tStartParallel);

% put parallel results together
FXYvec = zeros(size(Xvec));
for proc=1:numOfWorkers
    FXYvec(idxparallel{proc}) = FXYparallel{proc};
end

% back to matrix
FXY_Parallel = reshape(FXYvec,size(X));

% finalize parallel environment
delete(poolobj);

%% GIVE INFO
disp(['time sequential : ' num2str(timeSequential) 's'])
disp(['time parallel   : ' num2str(timeParallel) 's'])

% compute difference of results (should be zero)
disp(['diff = ' num2str(norm(FXY_Sequential - FXY_Parallel,1))])

%% PLOT RESULTS
% (just for fun)
%FXY = FXY_Sequential; % what to plot
FXY = FXY_Parallel; % what to plot

figure
hold on
mesh(X,Y,FXY)

hold off
