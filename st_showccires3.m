function st_showccires3(T,cx,ste)
   

    MASK=logical((cx{1}==1)*(cx{1}==2)');
    A=logical(full(sc_knngraph(ste.xy,10)));

    MASK1=MASK&A;
    MASK2=MASK&(~A);


    [x1,x2]=find(MASK1);
    a=ste.xy(x1,:);
    b=ste.xy(x2,:);
    ab=[mean([a(:,1) b(:,1)],2) mean([a(:,2) b(:,2)],2)];
   
%   group=ste.sce.struct_cell_clusterings.sc3;
%   idx1=3;
%   idx2=4;
% 
%      MASK=logical(group==idx1*(group==idx2)');
%      MASK1=MASK&A;
%      % isequal(MASK1,M)
%      MASK2=MASK&(~MASK1);    
     [i1,i2]=find(MASK1);    

methodid=1;
xx=1;
        [~,lidx]=ismember(T.lgene(xx),ste.sce.g);
        [~,ridx]=ismember(T.rgene(xx),ste.sce.g);
        x1=ste.sce.X(lidx,:)';       % ligand
        x2=ste.sce.X(ridx,:)';       % receptor
    
        switch methodid
            case 1
                X12=x1*x2';
                nv=X12(MASK1);
                vv=X12(MASK2);
            case 2
                x1=x1(i1)';
                x2=x2(i2)';
                nv=x1.*x2;
                x2=x2(randperm(size(x2,1)));
                vv=x1.*x2;
        end

        size(nv)
        size(vv)

        [~,~,~]=kstest2(nv,vv);
        figure
        cdfplot(nv); 
        hold on 
        cdfplot(vv)


