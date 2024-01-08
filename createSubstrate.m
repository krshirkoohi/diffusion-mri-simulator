function [I, disks] = createSubstrate(ra,N_ii,C,layers)

    pop = C.pop;
    sf = C.scale;
    spacing = C.spacing;
    alpha = C.alpha;
    beta = C.beta;

    % Gamma distribution
    gamma = gamrnd(alpha,beta,pop,1);
    radii = sort(gamma','descend')*sf;

    % Bubblebath parameters
    Bb.frameSize=[N_ii,N_ii];
    Bb.overlapType='relative';
    Bb.overlap=spacing-1;
    Bb.density=(pop/(N_ii));
    Bb.circSize=radii;
    Bb.nSizes=NaN;
    Bb.edgeType=3;
    Bb.maxCircsPerRad=1;
    Bb.supressWarning=true;
    [BbData, BbHandles, Bb] = bubblebath(Bb);

    % Output
    disks.centers = BbData(:,1:2);
    disks.radii = BbData(:,3);
    disks.count = length(disks.radii);
    difference = (disks.count-pop);
    if difference<0
        disks.surplus = 0;
        disks.missed = norm(difference);
    else
        disks.surplus = norm(difference);
        disks.missed = 0;
    end

    % Create mask from output
    I = uint8(createCirclesMask([N_ii N_ii],disks.centers+0.5*N_ii,disks.radii));
    if layers > 1
        % Create mask from output
        I = I + uint8(createCirclesMask([N_ii N_ii],disks.centers+0.5*N_ii,disks.radii*0.6));
        % 0.6 = g ratio factor
    end
end

