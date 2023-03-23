% Grouping detected line segments into power lines using graph theory

function W = clusterEdges(Q)
CC = bwconncomp(Q);
A = zeros(CC.NumObjects);
a=1; b=1;

for i=1:CC.NumObjects
    
    for j=1:CC.NumObjects
        [xi, yi] = ind2sub(CC.ImageSize, CC.PixelIdxList{i});
        [xj, yj] = ind2sub(CC.ImageSize, CC.PixelIdxList{j});
        if isempty(intersect(xi, xj)) && isempty(intersect(yi, yj))
            if xi(1)<xj(1)
                x=[xi; xj];
                y=[yi; yj];
            else
                x=[xj; xi];
                y=[yj; yi];
            end
            pi = polyfit(xi,yi,1);
            pj = polyfit(xj,yj,1);
            p = polyfit(x,y,1);
            f = polyval(p,x, 'r');
            err= mean(abs(y-f));
            [uCentre, vCentre, Ru, Rv, thetarad] = fitellipse(x,y);
                
            ecc = sqrt(1-(min(Ru, Rv)^2/max(Ru, Rv)^2));
            w = abs(a * exp(-b * err/ecc));
           
            if w>0.00001

                A(i,j)=1;
            end
       end

    end
end

A=A+diag(ones(1, length(A)));
S=sparse(tril(A));
[I, J, s] = find(S);
line = struct([]);
for i=1:length(s)
    pair = sort([I(i) J(i)]);
    n = 0;
    for j=1:numel(line)
        currentPair = line(j).segment;
        if any(ismember(pair, currentPair))
            n=j;
        end
    end
    if n == 0
        line(numel(line)+1).segment = pair;
    else
        line(n).segment = unique([line(n).segment pair]);
    end
end

% morphological properties
W = zeros(size(Q));
CC = bwconncomp(Q);
for i=1:numel(line)
    C=line(i).segment;
    LInd=[];
    
    for j=C
        LInd=vertcat(LInd, CC.PixelIdxList{j});
    end
    V(LInd)=0;
    [x, y] = ind2sub(CC.ImageSize, sort(LInd));
    p = polyfit(x,y,1);
    x0= linspace(min(x),max(x), 1000)';
    f = round(polyval(p,x0, 'r'));
    f(find(f<1))=1;
    f(find(isnan(f)))=1;
    f(find(f>CC.ImageSize(2)))=CC.ImageSize(2);
    x0=round(x0);
    [x0 f];
    W(sub2ind(CC.ImageSize,x0,f))=1;
    IA=setdiff(1:length(A),C);
    A = A(IA,IA);
end

W=bwmorph(Q, 'skel', inf);

bp=find(bwmorph(W, 'branchpoints'));
for i=1:numel(bp)
    [i,j] = ind2sub(size(Q), bp(i));
    if any(~[i-1:i+1,j-1:j+1])==0
        W(i-1:i+1,j-1:j+1)=0;
    end
end
W=W(1:size(Q,1),1:size(Q,2));