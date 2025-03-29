function [nodePairsData] = nodeManager(interf, coord, interfData)

    topNodes = interf(:,1:4);
    botNodes = interf(:,5:8);
    topNodes = topNodes(:);
    botNodes = botNodes(:);
    nodePairs = [topNodes, botNodes];
    nodePairs = unique(nodePairs, 'rows');
    nni = size(nodePairs, 1);
    nodePairsData = repmat(struct('ntop', 1, 'nbottom', 1, ...
        'normal', zeros(3,1), 'coord', zeros(3,1)), nni, 1);
    % Compute a representative normal for each node
    topFaces = interf(:,1:4);           % ni * 4
    for i = 1 : nni
        % extract the faces sharing node i
        nodePairsData(i).ntop = nodePairs(i,1);
        nodePairsData(i).nbottom = nodePairs(i,2);
        facesSharingNode = find(any(topFaces == nodePairs(i,1), 2));
        n = zeros(3,1);
        for f = 1 : size(facesSharingNode, 1)
            fi = facesSharingNode(f);
            n = n + interfData(fi).normal;
        end
        nodePairsData(i).normal = n./size(facesSharingNode, 1);
        nodePairsData(i).coord = coord(nodePairs(i,1), :)';
    end
