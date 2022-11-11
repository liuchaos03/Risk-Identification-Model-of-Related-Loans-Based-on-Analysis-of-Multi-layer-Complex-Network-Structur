% aver_C=zeros();
% aver_w=zeros();
% for k=1:6
% nodes=char(nodes_strcu{1,k});
% nodes=string(nodes);
% 
% tmp=string(char(link2(:,1)));
% tmk=string(char(link2(:,2)));
% link3=string();
% link3=[tmp,tmk];
% 
% GA=int8(zeros(length(nodes)));
% 
% for i=1:length(link3)
%     tmk=link3(i,1)==nodes;
%     tmp=link3(i,2)==nodes;
%     if sum(tmk)*sum(tmp)>0
%         GA(tmk,tmp)=1;
%         disp(i);
%     end
% end
% 
% disp(sum(sum(GA)));
% [~,C]=Clustering_Coefficient(GA);
% wbc = centrality(digraph(single(GA)),'betweenness');
% aver_C(k)=C;
% aver_w(k)=mean(wbc);
% end

% link4=cell({});
% link4(:,1)=raw_linkbank(:,3);
% link4(:,2)=raw_linkbank(:,6);
% money=raw_linkbank(:,8);
% 
% tmp=string(char(link4(:,1)));
% tmk=string(char(link4(:,2)));
% link4=string();
% link4=[tmp,tmk];
% 
% money=cell2mat(money);
% tmk=link4(:,1)==link4(:,2);
% link4(tmk,:)=[];
% link4=strrep(link4,' ','');
aver_P=zeros();
aver_R=zeros();
for k=1:6
nodes=char(nodes_strcu{1,k});
nodes=string(nodes);
GA=zeros(length(nodes));

for i=1:length(link4)
    tmk=link4(i,1)==nodes;
    tmp=link4(i,2)==nodes;
    key=tmk+tmp;
    if sum(key)>=2
        GA(tmk,tmp)=sum(money(i))/sum(money(link4(:,1)==nodes(tmk)));
%         disp(i);
    end
end

P=zeros();
tmp=zeros();
for i=1:length(nodes)
   tmp(i,1)=sum(money(link4(:,2)==nodes(i)));  %个人的贷款总量
end
for i=1:length(nodes)
   P(i,1)=min(tmp(i),(GA(i,:)*tmp));
end
for i=1:length(P)
    if tmp(i)==0 && P(i,1)~=0
        P(i,1)=1;
        continue;
    end
    if tmp(i)==0 && P(i,1)==0
        P(i,1)=0;
        continue;
    end
    P(i,1)=P(i,1)./tmp(i,1);
end
aver_P(k)=mean(P);

I=eye(length(GA));
one=ones(1,length(GA));
Rs=((I-0.5.*single(GA))^-1)*one';

key=find(max(Rs)==Rs);
GA(key(1),:)=[];
GA(:,key(1))=[];

I=eye(length(GA));
one=ones(1,length(GA));
Rs=((I-0.5.*single(GA))^-1)*one';
Rs=(one*Rs)^-1;
% aver_r(k)=r;
aver_R(k)=mean(Rs);

disp(sum(sum(k)));
end

