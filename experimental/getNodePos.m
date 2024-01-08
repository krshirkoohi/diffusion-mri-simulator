function nodePos = getNodePos(node,X,Q)
    Y = X';
    nodePos = [X(Q.ind(node)),Y(Q.ind(node))];
end