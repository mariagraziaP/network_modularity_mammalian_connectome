function [ciu,Aall,anull,A,ciall,gam_range] = get_MRCC_mammals(K,sam1,sam2,maxC,r1,r2,tau)
% 
%   INPUTS:
%       K:      mammals data -> [N*N*M]
%       sam1:   samples communities first loop
%       sam2:   samples communities second loop
%       maxC:   max num communities
%       r1:     initial gamma low
%       r2:     initial gamma high
%       tau:    tau parameter for consensus modules at the end 
%
%   OUTPUTS:
%       ciu:    consensus community partition
%       Aall:   agreement matrix
%       anull:  analytic null
%       A:      agree - null
%       ciall:  community assignments

if nargin < 7
   tau = 0.0 ; 
end

% performs a 'lite' version of MRCC, obtaining a CC (Aall) matrix and a set of
% consensus modules

N = size(K,1); % number of nodes
M = size(K,3); % number of mammals

% FIRST PASS - find gamma range where modules range from 2 to (just below) maxC;
% initial gamma range
gam_range = logspace(r1,r2,sam1);
G = length(gam_range);

ciall = zeros(N,G,M);
parfor g=1:length(gam_range)
    disp([num2str(g) ' round1'])
%     if ~mod(g,100) ; disp([ num2str(g) ' of ' num2str(length(gam_range))]) ; end
    for m=1:M
        [ci, ~] = community_louvain(K(:,:,m),gam_range(g),[]);
        ciall(:,g,m) = ci;
    end
end

% identify the limits of the range
for m=1:M
    ff = find((max(squeeze(ciall(:,:,m)))>1)&(max(squeeze(ciall(:,:,m)))<maxC));
    g1(m) = log10(gam_range(ff(1)));
    g2(m) = log10(gam_range(ff(end)));
%     disp([num2str(g1),' ',num2str(g2)])
end

g1 = g1(find(max(g1),1));
g2 = g2(find(min(g2),2));

gam_range = logspace(g1,g2,sam2);
G = length(gam_range);

clear ciall
ciall = zeros(N,M,G);
parfor g=1:G
    disp([num2str(g) ' round2'])
    for m=1:M
        %     if ~mod(g,1000) ; disp([ num2str(g) ' of ' num2str(G)]) ; end
        [ci,~] = community_louvain(K(:,:,m),gam_range(g),[]);
        ciall(:,m,g) = ci;
    end
end

% exclude spurious partitions that are outside of the range 2 to maxC
for m=1:M
    numm = max(squeeze(ciall(:,m,:)));
%     use = find((numm>1)&(numm<maxC));
    use(:,m) = (numm>1)&(numm<maxC);
end
commuse = find(sum(use,2)==M);
ciall = ciall(:,:,commuse);

for m=1:M

    Aall(:,:,m) = agreement(squeeze(ciall(:,m,:)))./length(commuse);

    % analytic null
    anull = 0;
    for cnum = 1:length(commuse)
        anull = anull+sum((histcounts(ciall(:,m,cnum))./N).*((histcounts(ciall(:,m,cnum))-1)./(N-1)));
    end
    anull = anull/length(commuse);
    A(:,:,m) = Aall(:,:,m) - anull;

    % consensus clustering
    ciu(:,m) = consensus_und(squeeze(A(:,:,m)),tau,100);

end



