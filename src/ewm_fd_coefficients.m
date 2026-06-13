function coeff = ewm_fd_coefficients(kind, order)
%EWM_FD_COEFFICIENTS 交错网格一阶导数差分系数。

if nargin < 2 || isempty(order)
    order = 5;
end

kind = lower(char(kind));

switch kind
    case {'standard', 'taylor', 'classic'}
        if order == 5
            coeff = [ ...
                1.211242675807458, ...
               -0.089721679689456, ...
                0.013842773437802, ...
               -0.001765659877271, ...
                0.000118679488701];
        else
            coeff = taylor_staggered_coefficients(order);
        end
    otherwise
        error('Unknown coefficient kind "%s".', kind);
end
end

function coeff = taylor_staggered_coefficients(order)
offset = (0:order-1) + 0.5;
A = zeros(order, order);
b = zeros(order, 1);

for row = 1:order
    power = 2 * row - 1;
    A(row, :) = 2 * offset .^ power;
    if row == 1
        b(row) = 1;
    end
end

coeff = (A \ b).';
end
