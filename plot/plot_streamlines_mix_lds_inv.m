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

% Plot dynamics
for c=1:n_comp
    weights(c,:) = ( mvnpdf(x', lambda.mu_xloc{c}', ...
                                         lambda.cov_xloc{c}) ...
                   .* lambda.pi(c) )' + realmin; % TODO: In case of numerical 
                                                 % problems, instead of realmin
                                                 % choose closest component 
                                                 % based on Mahalanobis dist
end
weights = weights ./ repmat(sum(weights,1)+n_comp*realmin, n_comp, 1); 

sum_A_inv = zeros(2,2,length(x));
for c=1:n_comp
    eig(inv(lambda.A_inv{c}))
    if (sum( eig(lambda.A_inv{c} + lambda.A_inv{c}') <= 0 ) > 0)
        disp('Not stable!')
    end
    sum_A_inv = sum_A_inv + repmat(reshape(weights(c,:), ...
                                [1 1 length(weights(c,:))]), 2,2,1)...
                                .*repmat(lambda.A_inv{c},1,1,length(x));
end
for i=1:length(x)
    x_dot(:,i) = -sum_A_inv(:,:,i)\(x(:,i) - lambda.x_attractor);
end
x_dyn_h = streamslice(x_tmp,y_tmp,reshape(x_dot(1,:),ny,nx), ...
                                reshape(x_dot(2,:),ny,nx),1,'method','cubic');
hold on;
% Plot attractor
x_attractor_h = plot(lambda.x_attractor(1), lambda.x_attractor(2), ...
                                        'o', 'LineWidth', 6,'MarkerSize', 12);
axis([ax.XLim ax.YLim]);
box on;
legend([x_dyn_h(1) x_attractor_h], 'xdot', 'attractor');

for c=1:n_comp
    p_handle = plot_ellipsoid(lambda.mu_xloc{c}, ...
                                        lambda.cov_xloc{c});
    set(p_handle, 'EdgeColor',[1 0.5 0], 'EdgeAlpha', 0.3, ...
                  'FaceColor',[1 0.5 0], 'FaceAlpha', 0.05);
end

legend([x_dyn_h(1) x_attractor_h p_handle], 'xdot', 'attractor', ...
                                                            'Local component');
end

