function plot_streamlines(A, b, limits)
d = size(A,1);
if d~=2
    disp('This function can only be used for 2D settings.')
    return
end

nx=100;
ny=100;

ax.XLim = limits(1:2);
ax.YLim = limits(3:4);
ax_x=linspace(ax.XLim(1),ax.XLim(2),nx); 
ax_y=linspace(ax.YLim(1),ax.YLim(2),ny); 
[x_tmp,y_tmp]=meshgrid(ax_x,ax_y); 
x=[x_tmp(:) y_tmp(:)]';

x_dot = A*x + repmat(b, [1 size(x,2)]); 
streamslice(x_tmp,y_tmp,reshape(x_dot(1,:),ny,nx),reshape(x_dot(2,:),ny,nx),1,'method','cubic')
axis([ax.XLim ax.YLim]);
box on;

end

