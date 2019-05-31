function [ mean_arr ] = make_mean( filename )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
full_arr=readfile4(filename);
[temp,ind]=sort(full_arr(:,1));
sorted=cat(2,full_arr(ind,1),full_arr(ind,2),full_arr(ind,3),full_arr(ind,4));

mean_arr=make_mean_for_arr(sorted,4);
% place=1;
% array_size=1;
% count=0;
% 
%     for i=2:length(x)      
%         if (x(i)~=x(i-1))
%             array_size=array_size+1;
%         end
%     end
% arr=zeros(array_size,4);
%     for i=1:length(x)      
%         if(i>1)
%             if (x(i)~=x(i-1))
%                 arr(place,:)=arr(place,:)/count;
%                 place=place+1;
%                 count=0;
%             end
%         end
%         arr(place,:)=arr(place,:)+full_arr(i,:);
%         count=count+1;
%     end
%     arr(place,:)=arr(place,:)/count;
     %mean_arr=arr;  
     
     filename=strcat(filename,'_mean');
     writefile(filename, mean_arr,4);
     
     
end

