function demo
    a = [1 3]';
    b = [1 1 2]';
    C = [1 2; 3 4; 5 6]';
    T = auction(a,b,-C);
    
    figure(1)
    createPlot(a,b,C,T,false)
    
    figure(2)
    createPlot(a,b,C,T,true)
end

function createPlot(a,b,C,T,plotSolution)
    sourceX = 0;
    sinkX  = 10;

    clf
    for i=1:length(a)
       text(sourceX, i, sprintf('source %d (cap %d)', i, a(i)), 'HorizontalAlignment','right');
       rectangle('Position', [sourceX-5, i-0.1,5, 0.2]);
       for j=1:length(b)
           if ~plotSolution
               line([sourceX sinkX], [i j])
               alpha = 0.7;
               text(alpha*(sinkX-sourceX), (1-alpha)*i + alpha*j, sprintf('cost %d', C(i,j)))
           else
                if T(i,j)>0
                    line([sourceX sinkX], [i j])
                    alpha = 0.7;
                    text(alpha*(sinkX-sourceX), (1-alpha)*i + alpha*j, sprintf('flow %d', T(i,j)))
                end
           end
       end
    end

    for i=1:length(b)
       text(sinkX, i, sprintf('sink %d (cap %d)', i, b(i)));
       rectangle('Position', [sinkX, i-0.1,4, 0.2]);
    end

    ylim([0.5, max(length(a), length(b))+0.5])
    xlim([sourceX-5 sinkX+5])
    axis off
end