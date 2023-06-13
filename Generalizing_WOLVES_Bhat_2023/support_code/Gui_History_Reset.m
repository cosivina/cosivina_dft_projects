%author: Ajaz Bhat ajaz.bhat@ubd.edu.bn

 % GUI initilaizations and looking history updation 
if (mode == 1 )    
        gui.init();
        delete(gui.visualizations{1}.axesHandle);
        gui.visualizations{1}.axesProperties{2} = [-t_max, 10];
        gui.visualizations{1}.plotProperties{1}{end} = 0:-1:-t_max+1;
        gui.visualizations{1}.plotProperties{2}{end} = 0:-1:-t_max+1;
        gui.visualizations{1}.plotProperties{3}{end} = 0:-1:-t_max+1;
        gui.visualizations{1}.plotProperties{4}{end} = 0:-1:-t_max+1;
        gui.visualizations{1}.plotProperties{5}{end} = 0:-1:-t_max+1;
        gui.visualizations{1}.init(gui.figureHandle);
end