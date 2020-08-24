function stats = rm_anova4(Y,S,F1,F2,F3,F4,FACTNAMES)
%
% function stats = rm_anova2(Y,S,F1,F2,FACTNAMES)
%
% Two-factor, within-subject repeated measures ANOVA.
% For designs with two within-subject factors.
%
% Parameters:
%    Y          dependent variable (numeric) in a column vector
%    S          grouping variable for SUBJECT
%    F1         grouping variable for factor #1
%    F2         grouping variable for factor #2
%    FACTNAMES  a cell array w/ two char arrays: {'factor1', 'factor2'}
%
%    Y should be a 1-d column vector with all of your data (numeric).
%    The grouping variables should also be 1-d numeric, each with same
%    length as Y. Each entry in each of the grouping vectors indicates the
%    level # (or subject #) of the corresponding entry in Y.
%
% Returns:
%    stats is a cell array with the usual ANOVA table:
%      Source / ss / df / ms / F / p
%
% Notes:
%    Program does not do any input validation, so it is up to you to make
%    sure that you have passed in the parameters in the correct form:
%
%       Y, S, F1, and F2 must be numeric vectors all of the same length.
%
%       There must be at least one value in Y for each possible combination
%       of S, F1, and F2 (i.e. there must be at least one measurement per
%       subject per condition).
%
%       If there is more than one measurement per subject X condition, then
%       the program will take the mean of those measurements.
%
% Aaron Schurger (2005.02.04)
%   Derived from Keppel & Wickens (2004) "Design and Analysis" ch. 18
%

%
% Revision history...
%
% 11 December 2009 (Aaron Schurger)
% 
% Fixed error under "bracket terms"
% was: expY = sum(Y.^2);
% now: expY = sum(sum(sum(MEANS.^2)));
%

stats = cell(4,5);

F1_lvls = unique(F1);
F2_lvls = unique(F2);
F3_lvls = unique(F3);
F4_lvls = unique(F4);
Subjs = unique(S);

a = length(F1_lvls); % # of levels in factor 1
b = length(F2_lvls); % # of levels in factor 2
c = length(F3_lvls);
d = length(F4_lvls);
n = length(Subjs); % # of subjects

INDS = cell(a,b,c,d,n); % this will hold arrays of indices
CELLS = cell(a,b,c,d,n); % this will hold the data for each subject X condition
MEANS = zeros(a,b,c,d,n); % this will hold the means for each subj X condition

% Calculate means for each subject X condition.
% Keep data in CELLS, because in future we may want to allow options for
% how to compute the means (e.g. leaving out outliers > 3stdev, etc...).
for i=1:a % F1
    for j=1:b % F2
        for l=1:c
            for m=1:d
                for k=1:n % Subjs
                    INDS{i,j,l,m,k} = find(F1==F1_lvls(i) & F2==F2_lvls(j) & S==Subjs(k)& F3==F3_lvls(l)  & F4==F4_lvls(m));
                    CELLS{i,j,l,m,k} = Y(INDS{i,j,l,m,k});
                    MEANS(i,j,l,m,k) = mean(CELLS{i,j,l,m,k});
                end
            end
        end
    end
end

% make tables (see table 18.1, p. 402)
ABCD = reshape(sum(MEANS,5),a,b,c,d); % across subjects
AS = reshape(sum(MEANS,[2 3 4]),a,n); % across factor 2
BS = reshape(sum(MEANS,[1 3 4]),b,n); % across factor 1
CS = reshape(sum(MEANS,[1 2 4]),c,n); % across factor 3
DS = reshape(sum(MEANS,[1 2 3]),d,n); % across factor 4

A = sum(ABCD,[2 3 4]); % sum across columns, so result is ax1 column vector
B = sum(ABCD,[1 3 4]); % sum across rows, so result is 1xb row vector
C = sum(ABCD,[1 2 4]);
D = sum(ABCD,[1 2 3]);

S = sum(AS,1); % sum across columns, so result is 1xs row vector
T = sum(sum(A)); % could sum either A or B or S, choice is arbitrary

% degrees of freedom
dfA = a-1;
dfB = b-1;
dfC = c-1;
dfD = d-1;

dfAB = (a-1)*(b-1);
dfAC = (a-1)*(c-1);
dfAD = (a-1)*(d-1);
dfBC = (b-1)*(c-1);
dfBD = (b-1)*(d-1);
dfCD = (c-1)*(d-1);

dfS = n-1;
dfAS = (a-1)*(n-1);
dfBS = (b-1)*(n-1);
dfCS = (c-1)*(n-1);
dfDS = (d-1)*(n-1);

dfABS = (a-1)*(b-1)*(n-1);
dfACS = (a-1)*(c-1)*(n-1);
dfADS = (a-1)*(d-1)*(n-1);
dfBCS = (b-1)*(c-1)*(n-1);
dfBDS = (b-1)*(d-1)*(n-1);
dfCDS = (c-1)*(d-1)*(n-1);

% bracket terms (expected value)
expA = sum(A.^2)./(b*c*d*n);
expB = sum(B.^2)./(a*c*d*n);
expC = sum(C.^2)./(a*b*d*n);
expD = sum(D.^2)./(a*b*c*n);

expAB = sum(sum(AB.^2))./[c*d*n];
expAC = sum(sum(AC.^2))./[b*d*n];
expAD = sum(sum(AD.^2))./[b*c*n];
expBC = sum(sum(BC.^2))./[a*d*n];
expBD = sum(sum(BD.^2))./[a*c*n];
expCD = sum(sum(CD.^2))./[a*b*n];

expS = sum(S.^2)./(a*b*c*d);

expAS = sum(sum(AS.^2))./[b*c*d];
expBS = sum(sum(BS.^2))./[a*c*d];
expCS = sum(sum(CS.^2))./[a*b*d];
expDS = sum(sum(DS.^2))./[a*b*c];

expY = sum(sum(sum(sum(sum(MEANS.^2))))); %sum(Y.^2);
expT = T^2 / (a*b*c*d*n);

% sums of squares
ssA = expA - expT;
ssB = expB - expT;
ssC = expC - expT;
ssD = expD - expT;

ssAB = expAB - expA - expB + expT;
ssAC = expAC - expA - expC + expT;
ssAD = expAD - expA - expD + expT;
ssBC = expBC - expB - expC + expT;
ssBD = expBD - expB - expD + expT;
ssCD = expCD - expC - expD + expT;

ssS = expS - expT;
ssAS = expAS - expA - expS + expT;
ssBS = expBS - expB - expS + expT;
ssCS = expCS - expC - expS + expT;
ssDS = expDS - expD - expS + expT;

ssABS = expY - expAB - expAS - expBS + expA + expB + expS - expT;
ssACS = expY - expAC - expAS - expCS + expA + expC + expS - expT;
ssADS = expY - expAD - expAS - expDS + expA + expD + expS - expT;
ssBCS = expY - expBC - expBS - expCS + expB + expC + expS - expT;
ssBDS = expY - expBD - expBS - expDS + expB + expD + expS - expT;
ssCDS = expY - expCD - expCS - expDS + expC + expD + expS - expT;

ssTot = expY - expT;

% mean squares
msA = ssA / dfA;
msB = ssB / dfB;
msC = ssC / dfC;
msD = ssD / dfD;

msAB = ssAB / dfAB;
msAC = ssAC / dfAC;
msAD = ssAD / dfAD;
msBC = ssBC / dfBC;
msBD = ssBD / dfBD;
msCD = ssCD / dfCD;

msS = ssS / dfS;

msAS = ssAS / dfAS;
msBS = ssBS / dfBS;
msCS = ssCS / dfCS;
msDS = ssDS / dfDS;

msABS = ssABS / dfABS;
msACS = ssACS / dfACS;
msADS = ssADS / dfADS;
msBCS = ssBCS / dfBCS;
msBDS = ssBDS / dfBDS;
msCDS = ssCDS / dfCDS;

% f statistic
fA = msA / msAS;
fB = msB / msBS;
fC = msC / msCS;
fD = msD / msDS;

fAB = msAB / msABS;
fAC = msAC / msACS;
fAD = msAD / msADS;
fBC = msBC / msBCS;
fBD = msBD / msBDS;
fCD = msCD / msCDS;

% p values
pA = 1-fcdf(fA,dfA,dfAS);
pB = 1-fcdf(fB,dfB,dfBS);
pC = 1-fcdf(fC,dfC,dfCS);
pD = 1-fcdf(fD,dfD,dfDS);

pAB = 1-fcdf(fAB,dfAB,dfABS);
pAC = 1-fcdf(fAC,dfAC,dfACS);
pAD = 1-fcdf(fAD,dfAD,dfADS);
pBC = 1-fcdf(fBC,dfBC,dfBCS);
pBD = 1-fcdf(fBD,dfBD,dfBDS);
pCD = 1-fcdf(fCD,dfCD,dfCDS);

% return values
stats = {'Source','SS','df','MS','F','p';...
         FACTNAMES{1}, ssA, dfA, msA, fA, pA;...
         FACTNAMES{2}, ssB, dfB, msB, fB, pB;...
         FACTNAMES{3}, ssC, dfC, msC, fC, pC;...
         FACTNAMES{4}, ssD, dfD, msD, fD, pD;...
         [FACTNAMES{1} ' x ' FACTNAMES{2}], ssAB, dfAB, msAB, fAB, pAB;...
         [FACTNAMES{1} ' x ' FACTNAMES{3}], ssAC, dfAC, msAC, fAC, pAC;...
         [FACTNAMES{1} ' x ' FACTNAMES{4}], ssAD, dfAD, msAD, fAD, pAD;...
         [FACTNAMES{2} ' x ' FACTNAMES{3}], ssBC, dfBC, msBC, fBC, pBC;...
         [FACTNAMES{2} ' x ' FACTNAMES{4}], ssBD, dfBD, msBD, fBD, pBD;...
         [FACTNAMES{3} ' x ' FACTNAMES{4}], ssCD, dfCD, msCD, fCD, pCD;...
         [FACTNAMES{1} ' x Subj'], ssAS, dfAS, msAS, [], [];...
         [FACTNAMES{2} ' x Subj'], ssBS, dfBS, msBS, [], [];...
         [FACTNAMES{1} ' x ' FACTNAMES{2} ' x Subj'], ssABS, dfABS, msABS, [], []};
 
 return