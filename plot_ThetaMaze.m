function plot_ThetaMaze(maze)
% Plotting a theta paze  based on maze parameters
rectangle('Position',[-maze.radius_out, -maze.radius_out, maze.radius_out*2, maze.radius_out*2],'Curvature',[1 1]), hold on
rectangle('Position',[-maze.radius_in, -maze.radius_in, maze.radius_in*2, maze.radius_in*2],'Curvature',[1 1])

plot([-maze.arm_half_width -maze.arm_half_width],[+maze.cross_radii -maze.cross_radii],'k')
plot([maze.arm_half_width maze.arm_half_width],[+maze.cross_radii -maze.cross_radii],'k')

if isfield(maze,'reward_points')
    for i = 1:length(maze.reward_points)
        [x,y] = pol2cart(maze.reward_points(i)/maze.radius_in-3*pi/2,[maze.radius_in,maze.radius_out]);
        plot(-x,y,'--ok')
    end
end
