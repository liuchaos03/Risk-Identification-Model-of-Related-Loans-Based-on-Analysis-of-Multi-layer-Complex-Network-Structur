%----有向图---------------

logtxt=cell({});
tic;
% %%---------读取数据--------------------------------
% [~,~,raw_nodes] = xlsread("F:\竞赛与项目\项目2020\论文3\data\1.贷款数据匿名化处理成品.xlsx");
% [~,~,raw_linkbank] = xlsread("F:\竞赛与项目\项目2020\论文3\data\2.担保数据匿名化处理成品.xlsx");
% [~,~,raw_linkman] = xlsread("F:\竞赛与项目\项目2020\论文3\data\3.高管、股东、对外投资数据匿名化处理成品.xlsx");
% raw_linkbank(1,:)=[];
% raw_nodes(1,:)=[];
% raw_linkman(1,:)=[];
%
% % %-------多层网络的合并----------------------------
% links=cell({{},{},{}});
% ti=1;
% for i=1:size(raw_linkbank,1)
%     guarantee=raw_linkbank(i,3);
%     borrowere=raw_linkbank(i,6);
%     if sum(isnan(guarantee{1,1}))~=0|| sum(isnan(borrowere{1,1}))~=0
%         logtxt(end+1)={["出现空值-","raw_linkbank:",num2str(i)]};
%         continue;
%     end
%     if  strcmp(guarantee{1, 1},'ActiveX VT_ERROR: ')||  strcmp(borrowere{1, 1},'ActiveX VT_ERROR: ')
%         logtxt(end+1)={["出现错误-","raw_linkbank:",num2str(i)]};
%         continue;
%     end
%     if i==1
%         links(ti,1)=guarantee;
%         links(ti,2)=borrowere;
%         links(ti,3)=raw_linkbank(i,1);
%         ti=ti+1;
%         continue;
%     end
%     if sum(ismember(links(:,1),guarantee{1,1}).*ismember(links(:,2),borrowere{1,1}))==1
%         k=find(ismember(links(:,1),guarantee{1,1}).*ismember(links(:,2),borrowere{1,1})==1);
%         tmp=links(k,3);
%         tmp(end+1)=raw_linkbank(i,1);
%         links(k,3)={tmp};
%         continue;
%     end
%     links(ti,1)=guarantee;
%     links(ti,2)=borrowere;
%     links(ti,3)=raw_linkbank(i,1);
%     ti=ti+1;
%     disp(i);
% end
%
% ti=size(links,1);
% for i=1:size(raw_linkman,1)
%     guarantee=raw_linkman(i,3);
%     borrowere=raw_linkman(i,6);
%     if sum(isnan(guarantee{1,1}))~=0|| sum(isnan(borrowere{1,1}))~=0
%         logtxt(end+1)={["出现空值-","raw_linkman:",num2str(i)]};
%         continue;
%     end
%     if  strcmp(guarantee{1, 1},'ActiveX VT_ERROR: ')||  strcmp(borrowere{1, 1},'ActiveX VT_ERROR: ')
%         logtxt(end+1)={["出现错误-","raw_linkman:",num2str(i)]};
%         continue;
%     end
%     if sum(ismember(links(:,1),guarantee{1,1}).*ismember(links(:,2),borrowere{1,1}))==1
%         k=find(ismember(links(:,1),guarantee{1,1}).*ismember(links(:,2),borrowere{1,1})==1);
%         tmp=links(k,3);
%         tmp(end+1)=raw_linkbank(i,7);
%         links(k,3)={tmp};
%         continue;
%     end
%     links(ti,1)=guarantee;
%     links(ti,2)=borrowere;
%     links(ti,3)=raw_linkbank(i,7);
%     ti=ti+1;
%     disp(i);
% end
%
% %%-------------自相关链路的去除--------------------------------
% link2=cell({{},{},{}});
% k=0;
% for i=1:length(links)
%     if  strcmp(links{i,1},links{i,2}) &&  (sum(ismember(links(:,1),links{i,1}))+sum(ismember(links(:,2),links{i,2})))<=2
%         continue;
%     else
%         k=k+1;
%         link2(k,1:3)=links(i,1:3);
%     end
%     disp(i);
% end
%
%
%%-------网络结构的计算----------------------------
network_Recipr=cell({});
network_Circle=cell({});
network_Outstar=cell({});
network_Intstar=cell({});

network_Chain=cell({});
network_Weak=cell({});

[G,nodelist,g]=links2G(link2);
[bins,binsizes] = conncomp(G);
nodelist=nodelist';
for i=1:length(binsizes)
    if binsizes(i)>1
        nodes=nodelist(bins==i);
        disp([num2str(i),': ',num2str(length(nodes))]);
        %----圈联网络的计算------------
        if binsizes(i)==2
            if ismember_cells(link2,nodes{1,1},nodes{1,2})
                network_Recipr(end+1)={nodes};
            end
        end
        
        if binsizes(i)==3
            key=0;
            if ismember_cell(link2,nodes{1,1},nodes{1,2}) &&  ismember_cell(link2,nodes{1,2},nodes{1,3})  &&  ismember_cell(link2,nodes{1,3},nodes{1,1}) ...
                    || ismember_cell(link2,nodes{1,2},nodes{1,1}) &&  ismember_cell(link2,nodes{1,3},nodes{1,2})  &&  ismember_cell(link2,nodes{1,1},nodes{1,3})
                network_Circle(end+1)={nodes};
                key=1;
            end
            
            if ismember_cells(link2,nodes{1,1},nodes{1,2}) && ismember_cells(link2,nodes{1,2},nodes{1,3}) ...
                    || ismember_cells(link2,nodes{1,1},nodes{1,2}) && ismember_cells(link2,nodes{1,3},nodes{1,1}) ...
                    || ismember_cells(link2,nodes{1,2},nodes{1,3}) && ismember_cells(link2,nodes{1,3},nodes{1,1})
                network_Recipr(end+1)={nodes};
                key=2;
            end
            
            if key==0
                [g,G]=Circle_G(link2,nodes);
                tmp=distances(g);tmp(tmp==inf)=0;tmp=max(max(tmp)); 
                if sum(sum(G))==length(nodes)-1 && tmp==length(nodes)-1
                    network_Chain(end+1)={nodes};      %链路网络
                end
            end
            
            if key~=2
                network_Weak(end+1)={nodes};   %弱链接网络
            end
            
        end
        
        if binsizes(i)>=4
            for j=1:length(nodes)
                tmk=ismember(link2(:,1),nodes{1,j});  %源网络
                if sum(tmk)>=3
                    tmp=ismember(link2(:,2),nodes{1,j});
                    if (sum(tmp)==0)
                        network_Outstar(end+1)={[nodes(1,j),link2(tmk,1)']};
                    end
                end
                
                tmk=ismember(link2(:,2),nodes{1,j});  %汇网络
                if sum(tmk)>=3
                    tmp=ismember(link2(:,1),nodes{1,j});
                    if (sum(tmp)==0)
                        network_Intstar(end+1)={[nodes(1,j),link2(tmk,2)']};
                    end
                end
            end
            
            
            Circle_tmp=cell({});
            [g,G]=Circle_G(link2,nodes);
            distance=distances(g);distance(distance==inf)=0;distance=max(max(distance)); 
            key=0;
            for j=1:length(nodes)
                for k=1:length(nodes)
                    if j~=k
                        if ismember_cells(link2,nodes{1,j},nodes{1,k})
                            network_Recipr(end+1)={[nodes(1,j),nodes(1,k)]};
                        end
                    end
                end
                tmp=dfsearch(g,j);
                for k=length(tmp):-1:3
                    if G(tmp(k),j)==1
                        tmk=tmp(1:k);
                        while length(tmk)>=2
                            for all=length(tmk):-1:2
                                if  G(tmk(all-1),tmk(all))~=1
                                    tmk(all-1)=[];
                                    break
                                end
                            end
                            if all==2
                                break
                            end
                        end
                        if length(tmk)>=3
                            Circle_tmp(end+1)={nodes(tmk)};
                            key=key+1;
                            %  break      识别多重圈
                        end
                    end
                end
            end
            %-------圈网络去重-------------------
            if key~=0
                if key==1
                    network_Circle(end+1)=Circle_tmp(1,1);
                else
                    tmp=zeros();
                    for j=length(Circle_tmp):2
                        for k=j+1:1
                            if length(Circle_tmp{1,j})==length(Circle_tmp{1,k})
                                if sum(ismember(Circle_tmp{1,j},Circle_tmp{1,k}))==length(Circle_tmp{1,k})
                                    tmp(end+1)=j;
                                    break;
                                end
                            end
                        end
                    end
                    if sum(tmp)~=0
                        Circle_tmp(tmp)=[];
                    end
                    for j=1:length(Circle_tmp)
                        network_Circle(end+1)=Circle_tmp(1,j);
                    end
                end
            else
                tmp=distances(g);tmp(tmp==inf)=0;tmp=max(max(tmp)); 
                if sum(sum(G))==length(nodes)-1 && tmp==length(nodes)-1
                    network_Chain(end+1)={nodes};      %链路网络
                end
            end
            network_Weak(end+1)={nodes};   %弱链接网络
        end
    end
end
toc;
%%-------------------计算不良率--------------------------------
nodes_Recipr=network_Recipr{1,1};
nodes_Circle=network_Circle{1,1};
nodes_Chain=network_Chain{1,1};
nodes_Outstar=network_Outstar{1,1};
nodes_Intstar=network_Intstar{1,1};
nodes_weak=network_Weak{1,1};

for i=2:length(network_Recipr)
    nodes_Recipr=[nodes_Recipr,network_Recipr{1,i}];
end
for i=2:length(network_Circle)
    nodes_Circle=[nodes_Circle,network_Circle{1,i}];
end
for i=2:length(network_Chain)
    nodes_Chain=[nodes_Chain,network_Chain{1,i}];
end
for i=2:length(network_Outstar)
    nodes_Outstar=[nodes_Outstar,network_Outstar{1,i}];
end
for i=2:length(network_Intstar)
    nodes_Intstar=[nodes_Intstar,network_Intstar{1,i}];
end
for i=2:length(network_Weak)
    nodes_weak=[nodes_weak,network_Weak{1,i}];
end
nodes_Recipr=unique(nodes_Recipr);
nodes_Circle=unique(nodes_Circle);
nodes_Chain=unique(nodes_Chain);
nodes_Outstar=unique(nodes_Outstar);
nodes_Intstar=unique(nodes_Intstar);
nodes_weak=unique(nodes_weak);
%%--合并--
nodes_strcu={nodes_Recipr,nodes_Circle,nodes_Chain,nodes_Outstar,nodes_Intstar,nodes_weak};
%%--计算--
Output_bad=zeros();
for i=1:length(nodes_strcu)
    tmp=nodes_strcu{1,i};
    key=0;
    bad=0;
    all=0;
    for j=1:length(tmp)
        tmk=ismember(raw_nodes(:,3),tmp{1,j});
        if sum(tmk)==0
            key=key+1;
        else
            nodes=raw_nodes(tmk,:);
            all=all+sum(cell2mat(nodes(:,5)));
            tmk=ismember(nodes(:,6),'损失  ');
            bad=bad+sum(cell2mat(nodes(tmk,5)));
            tmk=ismember(nodes(:,6),'次级  ');
            bad=bad+sum(cell2mat(nodes(tmk,5)));
            tmk=ismember(nodes(:,6),'可疑  ');
            bad=bad+sum(cell2mat(nodes(tmk,5)));
        end
    end
    Output_bad(i)=bad/all;
end
%%-----计算随机-----
net_rand=cell({});
for i=1:length(nodes_strcu)
    nodes=length(nodes_strcu{1,i});
    all=zeros();
    for j=1:10000        %进行一万次等规模随机抽样
        tmp= randperm(length(raw_nodes),nodes);
        bad=0;
        tmk=ismember(raw_nodes(tmp,6),'损失  ');
        bad=bad+sum(cell2mat(raw_nodes(tmp(tmk),5)));
        tmk=ismember(raw_nodes(tmp,6),'次级  ');
        bad=bad+sum(cell2mat(raw_nodes(tmp(tmk),5)));
        tmk=ismember(raw_nodes(tmp,6),'可疑  ');
        bad=bad+sum(cell2mat(raw_nodes(tmp(tmk),5)));
        all(j)=bad/sum(cell2mat(raw_nodes(tmp,5)));
    end
    net_rand(1,i)={all};
end

function [g,G]=Circle_G(links,nodes)
G=zeros(length(nodes));
for i=1:length(nodes)
    sets=getfriends_s(links,nodes{1,i});
    key=ismember(nodes,sets);
    G(i,key)=1;
end
g=digraph(G);
end

function sets=getfriends_s(links,node)
% sets=[];
tmp=ismember(links(:,1),node);
sets=links(tmp,2);
% set1=links(tmp,2);
% 
% tmp=ismember(links(:,2),node);
% set2=links(tmp,1);
% 
% if ~isempty(set1)
%     sets=set1;
% end
% if  ~isempty(set2)
%     sets=[sets;set2];
% end
% if ~isempty(sets)
%     tmp=ismember(sets,node);
%     sets(tmp)=[];  %去除自己
% end
end

function sets=getfriends(links,node)
sets=[];
tmp=ismember(links(:,1),node);
% sets=links(tmp,2);
set1=links(tmp,2);

tmp=ismember(links(:,2),node);
set2=links(tmp,1);

if ~isempty(set1)
    sets=set1;
end
if  ~isempty(set2)
    sets=[sets;set2];
end
if ~isempty(sets)
    tmp=ismember(sets,node);
    sets(tmp)=[];  %去除自己
end
end



function [key,list,all]=findfriends(links,node,target,finded,list)

sets=getfriends(links,node);
if isempty(sets)
    key=0;
    list=[];
    all=[];
    return;
end

tmp=ismember(sets,finded);
sets(tmp)=[];
if isempty(sets)
    return;
end

if isempty(finded)
    sets(strcmp(sets,target{1,1}))=[];  %第一次不算
    finded=node;
end

if ismember(sets,target{1,1})
    list=[list,target];
    key=1;
    all(end+1)=list;
    return;
else
    key=0;
    finded=[finded,sets];
    finded=unique(finded);
    for i=1:length(sets)
        list=[list,node];
        [key,list,all]=findfriends(links,sets{i,1},finded,list);
    end
end

end

function [G,nodelist,g]=links2G(links)
nodelist=unique([links(:,1);links(:,2)]);
nodelist=unique(nodelist);
[~,tmpi]=ismember(links(:,1),nodelist);
[~,tmpj]=ismember(links(:,2),nodelist);
one=ones(length(links),1);
g=sparse([tmpi;tmpj],[tmpj,tmpi],[one;one],length(nodelist),length(nodelist));
nodenames=string(num2cell(1:length(nodelist)));
G=graph(g,nodenames);
end


function re=ismember_cell(links,node1,node2)
re=1==0;
tmp=ismember(links(:,1),node1);
tmk=ismember(links(:,2),node2);

if sum(tmp.*tmk)>=1
    re=1==1;
end

end


function re=ismember_cells(links,node1,node2)
re=1==0;
tmp=ismember(links(:,1),node1);
tmk=ismember(links(:,2),node2);

if sum(tmp.*tmk)>=1
    tmp=ismember(links(:,2),node1);
    tmk=ismember(links(:,1),node2);
    if sum(tmp.*tmk)>=1
        re=1==1;
    end
end

end
