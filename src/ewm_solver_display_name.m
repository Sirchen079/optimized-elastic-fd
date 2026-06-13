function name = ewm_solver_display_name(tag)
%EWM_SOLVER_DISPLAY_NAME 求解器的中文显示名称。

switch tag
    case 'regular_noabsorb'
        name = '常规网格（无吸收边界）';
    case 'staggered_noabsorb_standard'
        name = '交错网格（标准系数，无吸收边界）';
    case 'staggered_pml_standard'
        name = '交错网格（标准系数，PML）';
    case 'staggered_pml_minimax'
        name = '交错网格（基于最大范数目标函数的优化系数，PML）';
    case 'exp1_regular_noabsorb'
        name = '实验1：常规网格（无吸收边界）';
    case 'exp1_staggered_noabsorb_standard'
        name = '实验1：交错网格（标准系数，无吸收边界）';
    case 'exp2_staggered_noabsorb_standard'
        name = '实验2：交错网格（标准系数，无吸收边界）';
    case 'exp2_staggered_pml_standard'
        name = '实验2：交错网格（标准系数，PML）';
    case 'exp3_staggered_pml_standard'
        name = '实验3：交错网格（标准系数，PML）';
    case 'exp3_staggered_pml_minimax'
        name = '实验3：交错网格（基于最大范数目标函数的优化系数，PML）';
    otherwise
        name = tag;
end
end
