function h = plot_filled_arrow(x2, y2, lw, arrow_color)
% plot_filled_arrow Draw a 2D arrow with a filled triangular head.
%
% Inputs: 
%   x2: 1x2 vector [x_start x_end]
%   y2: 1x2 vector [y_start y_end]
%   lw: line width of the arrow shaft
%   arrow_color: MATLAB color specification (e.g. 'k' or [0 0 0])
%
% Output:
%   h: graphics handles [shaft_handle head_handle]

    % Ensure inputs are row vectors
    x2 = reshape(x2, 1, numel(x2));
    y2 = reshape(y2, 1, numel(y2));

    % Direction vector from start to end
    dx = x2(2) - x2(1);
    dy = y2(2) - y2(1);
    L = hypot(dx, dy);

    if L == 0
        warning('plot_filled_arrow: Zero length arrow requested.');
        h = gobjects(1, 2);
        return
    end

    ux = dx / L;
    uy = dy / L;

    % Arrowhead size proportional only to line width

    base_head_length = 0.69 * lw;
    base_head_width = 0.23 * lw;

    % For very short arrows, cap head length to avoid degeneracy
    max_head_length = 0.3 * L; max_head_width = 0.3 * L;
    head_length = min(base_head_length, max_head_length);
    head_width = min(base_head_width, max_head_width);

    % Tip point
    x_tip = x2(2);
    y_tip = y2(2);

    % Base (center) of the head
    x_base = x_tip - head_length * ux;
    y_base = y_tip - head_length * uy;

    % Perpendicular unit vector for head width
    nx = -uy;
    ny = ux;

    % Triangle vertices for the head
    x_head = [x_tip, x_base + 0.5 * head_width * nx, x_base - 0.5 * head_width * nx];
    y_head = [y_tip, y_base + 0.5 * head_width * ny, y_base - 0.5 * head_width * ny];

    % Plot shaft from start to base of head
    hold_state = ishold;
    hold on

    h(1) = plot([x2(1), x_base], [y2(1), y_base], 'LineWidth', lw, 'Color', arrow_color);

    % Plot filled arrowhead
    h(2) = patch(x_head, y_head, arrow_color, 'EdgeColor', arrow_color);

    if ~hold_state
        hold off
    end
end
