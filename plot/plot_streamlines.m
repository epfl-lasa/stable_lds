function [ output_args ] = plot_streamlines(A, b)
d = size(Sigma,1);
if d/2~=2
    disp('This function can only be used for 2D models!')
    return
end

quality='low';
if ~isempty(varargin)
	quality=varargin{1};
end

if strcmpi(quality,'high')
    nx=600;
    ny=600;
elseif strcmpi(quality,'medium')
    nx=400;
    ny=400;
else
    nx=200;
    ny=200;
end

ax.XLim = D(1:2);
ax.YLim = D(3:4);
ax_x=linspace(ax.XLim(1),ax.XLim(2),nx); %computing the mesh points along each axis
ax_y=linspace(ax.YLim(1),ax.YLim(2),ny); %computing the mesh points along each axis
[x_tmp y_tmp]=meshgrid(ax_x,ax_y); %meshing the input domain
x=[x_tmp(:) y_tmp(:)]';
z=zeros(1,nx*ny);

x_dot = repmat(A, size(x,1))*x + repmat(b, size(x,1)); 
streamslice(x_tmp,y_tmp,reshape(xd(1,:),ny,nx),reshape(xd(2,:),ny,nx),1,'method','cubic')
axis([ax.XLim ax.YLim]);box on

end

