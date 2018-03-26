% rng(3,'twister');
rownumber = 100000;
FinalMatrix = zeros(rownumber,12);
for overiter = 1:rownumber
alpha = rand(); %0.3;
beta  = rand(); %0.8;
gamma = rand(); %0.8;
delta = 1 + rand(); %1.5;
k = randi([10 20]); %10;
TargetSize = randi([10 4000]); %190; % Number of users that will be targeted
NoCastSel = randi([2 15]); %10; % Number of actors that will be cast in the final movie
NoDirSel = 1; % Number of directors that will be employed in the final movie
NoGenSel = randi([2 4]); %3; % Number of genres that will define the movie

MovieIDNameYearCreditsGenreString = ProjectMovieAll('ProjectMovieCredits.xlsx','Sheet1','A1:L240');

Cast = MovieIDNameYearCreditsGenreString(:,5:9);
CastCategorical = categorical(Cast);
CCT = CastCategorical';
CastIDs = vec2mat(grp2idx(CCT(:)),5);
MovieCastAdj = zeros(size(CastIDs,1),max(max(CastIDs)));
for ii = 1:size(CastIDs,1)
    MovieCastAdj(ii,CastIDs(ii,1:4)) = 1;
end
for ii = 1:size(CastIDs,1)
    if ~isnan(CastIDs(ii,5))
        MovieCastAdj(ii,CastIDs(ii,5)) = 1;
    end
end

Genre = MovieIDNameYearCreditsGenreString(:,10:12);
GenreCategorical = categorical(Genre);
GCT = GenreCategorical';
GenreIDs = vec2mat(grp2idx(GCT(:)),3);
MovieGenreAdj = zeros(size(GenreIDs,1),max(max(GenreIDs)));
for ii = 1:size(GenreIDs,1)
    for jj = 1:size(GenreIDs,2)
        if ~isnan(GenreIDs(ii,jj))
            MovieGenreAdj(ii,GenreIDs(ii,jj)) = 1;
        end
    end
end

DirectorCategorical = categorical(MovieIDNameYearCreditsGenreString(:,4));
DirectorIDs = grp2idx(DirectorCategorical);
MovieDirectorAdj = zeros(size(DirectorIDs,1),max(DirectorIDs));
for ii = 1:size(DirectorIDs,1)
    MovieDirectorAdj(ii,DirectorIDs(ii)) = 1;
end

FeatureMatrix = horzcat(MovieCastAdj,MovieDirectorAdj,MovieGenreAdj);

UserMovieRating = ProjectUserMovieRating('Projectfilteredallold.xls','users_filtered',1,16121);
UserIDCategorical = UserMovieRating.UserID;
UserIDNum = grp2idx(UserIDCategorical);
UserNumIDMovieRating = UserMovieRating;
UserNumIDMovieRating.UserID = UserIDNum;
UMRMatrix = table2array(UserNumIDMovieRating);
UMRMatrix(:,3) = 3+2*UMRMatrix(:,3);
[C,ia,ic] = unique(UMRMatrix(:,1:2),'rows');
meanval = accumarray(ic,UMRMatrix(:,3),[],@mean);
UMRMatrixUnique = [C, meanval];
RatingsMatrix = zeros(max(UMRMatrixUnique(:,2)),max(UMRMatrixUnique(:,1))); %movies-by-users!
for ii = 1:size(UMRMatrixUnique,1)
    RatingsMatrix(UMRMatrixUnique(ii,2),UMRMatrixUnique(ii,1)) = UMRMatrixUnique(ii,3);
end
OriginalRatingsMatrix = RatingsMatrix;
AvgRatingsMatrix = sum(RatingsMatrix)./sum(RatingsMatrix~=0);
for ci = 1:size(RatingsMatrix,2)
    RatingsMatrix(RatingsMatrix(:,ci)~=0,ci) = RatingsMatrix(RatingsMatrix(:,ci)~=0,ci)-AvgRatingsMatrix(ci);
end
% RatingsMatrix = sparse(RatingsMatrix);
PosRatingsMatrix = RatingsMatrix;
PosRatingsMatrix(PosRatingsMatrix<0) = 0;
NegRatingsMatrix = -RatingsMatrix;
NegRatingsMatrix(NegRatingsMatrix<0) = 0;
PosAdjRatMatrix = PosRatingsMatrix;
PosAdjRatMatrix(PosAdjRatMatrix>0) = 2.^(delta*PosAdjRatMatrix(PosAdjRatMatrix>0));
NegAdjRatMatrix = NegRatingsMatrix;
NegAdjRatMatrix(NegAdjRatMatrix>0) = 2.^(delta*NegAdjRatMatrix(NegAdjRatMatrix>0));
AdjRatingsMatrix = PosAdjRatMatrix + NegAdjRatMatrix;

FeatureProb_FM = FeatureMatrix./sum(FeatureMatrix); %movies-by-features
FeatureProb_FM = FeatureProb_FM'; %feature-by-movies
FeatureProb_MF = FeatureMatrix./sum(FeatureMatrix,2); %movies-by-features
DislikeProb_MU = NegAdjRatMatrix./sum(NegAdjRatMatrix,2); %movies-by-users
DislikeProb_MU(isnan(DislikeProb_MU)) = 0;
DislikeProb_UM = NegAdjRatMatrix./sum(NegAdjRatMatrix); %movies-by-users
DislikeProb_UM(isnan(DislikeProb_UM)) = 0;
DislikeProb_UM = DislikeProb_UM'; %users-by-movies
LikeProb_MU = PosAdjRatMatrix./sum(PosAdjRatMatrix,2); %movies-by-users
LikeProb_MU(isnan(LikeProb_MU)) = 0;
LikeProb_UM = PosAdjRatMatrix./sum(PosAdjRatMatrix); %movies-by-users
LikeProb_UM(isnan(LikeProb_UM)) = 0;
LikeProb_UM = LikeProb_UM'; %users-by-movies

DislikeWalk_FMU = FeatureProb_FM*DislikeProb_MU;
LikeWalk_FMU = FeatureProb_FM*LikeProb_MU;
DislikeWalk_FMUMU = FeatureProb_FM*DislikeProb_MU*DislikeProb_UM*DislikeProb_MU;
LikeWalk_FMUMU = FeatureProb_FM*LikeProb_MU*LikeProb_UM*LikeProb_MU;
DislikeWalk_FMFMU = FeatureProb_FM*FeatureProb_MF*FeatureProb_FM*DislikeProb_MU;
LikeWalk_FMFMU = FeatureProb_FM*FeatureProb_MF*FeatureProb_FM*LikeProb_MU;

FMU = LikeWalk_FMU-DislikeWalk_FMU; %2-step
FMUMU = LikeWalk_FMUMU-DislikeWalk_FMUMU; %User-based 4-step
FMFMU = LikeWalk_FMFMU-DislikeWalk_FMFMU; %Feature-based 4-step

PreferenceMatrix = alpha*FMU + beta*FMUMU + gamma*FMFMU; %feature-by-users

Uprime = randi([1 size(RatingsMatrix,2)],1,TargetSize);% Number of Users = size(RatingsMatrix,2)
SelPref = PreferenceMatrix(:,Uprime);
SelPref_F = sum(SelPref,2);
[~, SelCastIdx] = sort(SelPref_F(1:size(MovieCastAdj,2)),'descend'); % Number of Cast = size(MovieCastAdj,2)
[~, SelDirIdx] = sort(SelPref_F(1+size(MovieCastAdj,2):size(MovieDirectorAdj,2)+size(MovieCastAdj,2)),'descend'); % Number of Directors = size(MovieDirectorAdj,2)
[~, SelGenIdx] = sort(SelPref_F(1+size(MovieDirectorAdj,2)+size(MovieCastAdj,2):size(MovieGenreAdj,2)+size(MovieDirectorAdj,2)+size(MovieCastAdj,2)),'descend'); % Number of Genres = size(MovieGenreAdj,2)
SelDirIdx = SelDirIdx + size(MovieCastAdj,2);
SelGenIdx = SelGenIdx + size(MovieCastAdj,2) + size(MovieDirectorAdj,2);
SelMovieFeatures = zeros(1,size(FeatureMatrix,2));
SelMovieFeatures([SelCastIdx(1:NoCastSel);SelDirIdx(1:NoDirSel);SelGenIdx(1:NoGenSel)]) = 1;

UprimeRatingsMatrix = OriginalRatingsMatrix(:,Uprime);

v = sum(UprimeRatingsMatrix ~= 0,2)'; % Vector (number of movies)
Popular = v*FeatureMatrix; % POPULAR Baseline (Should narrow down the user selection to U')

r = sum(UprimeRatingsMatrix,2)./sum(UprimeRatingsMatrix~=0,2); % Vector (number of movies)
r(isnan(r)) = 0; r = r';
Ft = normc(FeatureMatrix);
Top = r*Ft; % TOP Baseline (Should narrow down the user selection to U')

[NeighborMovies_New,~] = knnsearch(FeatureMatrix,SelMovieFeatures,'k',k,'distance','cosine');
MovieRating_New = (1/length(NeighborMovies_New))*sum(r(NeighborMovies_New));
[NeighborMovies_Top,~] = knnsearch(FeatureMatrix,Top,'k',k,'distance','cosine');
MovieRating_Top = (1/length(NeighborMovies_Top))*sum(r(NeighborMovies_Top));
[NeighborMovies_Pop,~] = knnsearch(FeatureMatrix,Popular,'k',k,'distance','cosine');
MovieRating_Popular = (1/length(NeighborMovies_Pop))*sum(r(NeighborMovies_Pop));

%alpha - beta - gamma - delta - k - TargetSize - NoCastSel - NoDirSel -
%NoGenSel - MovieRating_New - MovieRating_Top - MovieRating_Popular
FinalMatrix(overiter,:) = [alpha, beta, gamma, delta, k, TargetSize, NoCastSel, NoDirSel, NoGenSel, MovieRating_New, MovieRating_Top, MovieRating_Popular];
end