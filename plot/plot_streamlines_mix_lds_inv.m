function plot_streamlines_mix_lds_inv(lambda, limits)
n_comp = length(lambda.pi);
d = size(lambda.A_inv{1},1);
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

weights = zeros(n_comp, size(x,2));
for c=1:n_comp
    weights(c,:) = ( mvnpdf(x', lambda.mu_xloc{c}', ...
                                         lambda.cov_xloc{c}) ...
                   .* lambda.pi(c) )';
end
weights = weights ./ repmat(sum(weights,1)+1e-30, n_comp, 1); 

sum_A_inv = zeros(2,2,length(x));
for c=1:n_comp
    sum_A_inv = sum_A_inv + repmat(reshape(weights(c,:), ...
                                [1 1 length(weights(c,:))]), 2,2,1)...
                                .*repmat(lambda.A_inv{c},1,1,length(x));
end

for i=1:length(x)
    x_dot(:,i) = -sum_A_inv(:,:,i)\(x(:,i) - lambda.x_attractor); 
end
streamslice(x_tmp,y_tmp,reshape(x_dot(1,:),ny,nx),reshape(x_dot(2,:),ny,nx),1,'method','cubic')
hold on;
plot(lambda.x_attractor(1), lambda.x_attractor(2), 'o', 'LineWidth', 6,'MarkerSize', 12);
axis([ax.XLim ax.YLim]);
box on;

for c=1:n_comp
    p_handle = plot_ellipsoid(lambda.mu_xloc{c}, ...
                                        lambda.cov_xloc{c});
    set(p_handle, 'EdgeColor',[1 0.5 0], 'EdgeAlpha', 0.3, ...
                  'FaceColor',[1 0.5 0], 'FaceAlpha', 0.05);
end

end

