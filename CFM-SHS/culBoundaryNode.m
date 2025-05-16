function boundaryNode = culBoundaryNode(adj,groundtruth)
n=numnodes(adj);
boundaryNode=[];
for i=1:n
    for j=1:n
        if (adj(i,j)==1&&(groundtruth(i)~=groundtruth(j)))
            boundaryNode=[boundaryNode,i];
            break;
        end
    end
end