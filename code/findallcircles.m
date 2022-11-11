function [circles, n_circle]=findallcircles(G,startid)
%Input:undirected graph G , startid
%Output:
n = size(G.Nodes,1);
A = adjacency(G);
strstart = num2str(startid);
gid_start = findnode(G,strstart);
id_v = dfsearch(G,gid_start);
str_v = dfsearch(G,strstart);
circles = cell(n);
n_circle = 0;
for j = 1:n
    circles{j,1} = strstart;
end
lencircle = 1;
for i = 1:length(id_v)-1
    if A(id_v(i),id_v(i+1))==1
        circles(n_circle+1,lencircle+1) = str_v(i+1);
        lencircle = lencircle+1;
        if A(gid_start,id_v(i+1))==1 && lencircle>2
            n_circle = n_circle+1;
            circles(n_circle+1,:) = circles(n_circle,:);
%             lencircle = 1;
        end
        
    elseif A(gid_start,id_v(i+1))==1
        circles(n_circle+1,:) = [];
        circles{n_circle+1,1} = strstart;
        circles(n_circle+1,2) = str_v(i+1);
        lencircle = 2;
    else
        1;
%         error('circle finding unexpected error!')
    end
end