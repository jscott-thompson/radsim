function [Thresh,NumCells]=getCFARLevel(M,Nt,Ng)

[Nr,Nc]=size(M);
Thresh=zeros(Nr,Nc);
NumCells=zeros(Nr,Nc);

for i=1:Nr
    for j=1:Nc
        
        Mtest=0;    Mguard=0;
        
        if i>Nt+Ng & j>Nt+Ng & i<Nr-Nt-Ng & j<Nc-Nt-Ng
            %   center case
            Mtest=M(i-Nt-Ng:i+Nt+Ng,j-Nt-Ng:j+Nt+Ng);
            Mguard=M(i-Ng:i+Ng,j-Ng:j+Ng);
        elseif i>Nt+Ng & i<Nr-Nt-Ng & j<Nc-Nt-Ng
            %   left edge case
            Mtest=M(i-Nt-Ng:i+Nt+Ng,1:j+Nt+Ng);
            Mguard=M(i-Ng:i+Ng,1:j+Ng); 
        elseif i>Nt+Ng & j>Nt+Ng & i<Nr-Nt-Ng
            %   right edge case
            Mtest=M(i-Nt-Ng:i+Nt+Ng,j-Nt-Ng:end);
            Mguard=M(i-Ng:i+Ng,j-Ng:end);
        elseif j>Nt+Ng & i<Nr-Nt-Ng & j<Nc-Nt-Ng
            %   top edge case
            Mtest=M(1:i+Nt+Ng,j-Nt-Ng:j+Nt+Ng);
            Mguard=M(1:i+Ng,j-Ng:j+Ng);
        elseif i>Nt+Ng & j>Nt+Ng & j<Nc-Nt-Ng
            %   bottom edge case
            Mtest=M(i-Nt-Ng:end,j-Nt-Ng:j+Nt+Ng);
            Mguard=M(i-Ng:end,j-Ng:j+Ng);            
        elseif i<Nr-Nt-Ng & j<Nc-Nt-Ng
            %   top left corner case
            Mtest=M(1:i+Nt+Ng,1:j+Nt+Ng);
            Mguard=M(1:i+Ng,1:j+Ng);            
        elseif j>Nt+Ng & i<Nr-Nt-Ng
            %   top right corner case
            Mtest=M(1:i+Nt+Ng,j-Nt-Ng:end);
            Mguard=M(1:i+Ng,j-Ng:end);
        elseif i>Nt+Ng & j<Nc-Nt-Ng
            %   bottom left corner case
            Mtest=M(i-Nt-Ng:end,1:j+Nt+Ng);
            Mguard=M(i-Ng:end,1:j+Ng);            
        elseif i>Nt+Ng & j<Nc-Nt-Ng
            %   bottom right corner case
            Mtest=M(i-Nt-Ng:end,j-Nt-Ng:end);
            Mguard=M(i-Ng:end,j-Ng:end);  
        end
            
        sumt=sum(Mtest(:)); nt=numel(Mtest);
        sumg=sum(Mguard(:)); ng=numel(Mguard); 
        
        if sumt<sumg
%             pause
        end
        
        Thresh(i,j)=(sumt-sumg)/(nt-ng);
        NumCells(i,j)=nt-ng;
    end
end

            