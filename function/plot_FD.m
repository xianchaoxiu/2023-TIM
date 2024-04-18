function plot_FD(statistic,threshold,line_width,font_size,det)
    N_fault = length(statistic);
    switch det
        case 1
            semilogy(statistic,'linewidth',line_width);
        case 2
            plot(statistic,'linewidth',line_width);
    end
    hold on
    plot([1 N_fault],[threshold threshold],'r--','linewidth',line_width);
    legend('Test statistic','Control limit','Location','NorthEast');
    set(gca,'FontSize',font_size);
end