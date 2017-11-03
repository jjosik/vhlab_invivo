function [ refit_mvg,p_array,opt_array,x_track_range,y_track_range,mu,sigma,e_array ] = ...
    mvn_density_singleellipse( kernel,s_kernel,...
    x,y,xmean,ymean,peak_val,watchfit )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

global scale;

%x_track_range = 10.^[-2:1:2];
x_track_range = [10:5:25];
y_track_range = 10.^[2:1:5];
%y_track_range = [0.02:0.02:0.2];
%x_track = x_track_range(randperm(length(x_track_range)));
p_array = cell(length(x_track_range),length(y_track_range));
opt_array = zeros(length(x_track_range),length(y_track_range));
e_array = cell(length(x_track_range),length(y_track_range));
for u = 1:length(y_track_range),
    for t = 1:length(x_track_range),
        cov0 = [x_track_range(t) 5; 5 y_track_range(u)];
        ch0 = chol(inv(cov0),'lower');
        params0 = [0,xmean,ymean,peak_val,ch0(1),ch0(3),ch0(4)];
        upper = [Inf;max(x(:));max(y(:));255;Inf;Inf;Inf];
        lower = [0;min(x(:));min(y(:));0;1e-15;-Inf;1e-15];
        opt = optimset('lsqnonlin');
        opt.Display = 'iter';
        opt.TolFun = 1e-8;
        %opt.Algorithm = 'levenberg-marquardt';
        presort = 1;
        [params,resnorm,resid,exitflag,output] = ...
            lsqnonlin(@(P) mvn2d_objfun(P,x,y,kernel,s_kernel,presort),params0,lower,upper,opt);
        p_array{t,u} = params;
        opt_array(t,u) = output.firstorderopt;
        cholcov = [params(5) 0; params(6) params(7)];
        covm = cholcov*cholcov';
        mu0 = [params0(2) params0(3)];
        mu = [params(2) params(3)];
        sigma = inv(covm);
        sigma0 = cov0;
        if watchfit == 1,
            figure(100);
            subplot(2,3,1); imagesc(kernel); colorbar; title('Original kernel');   %1. original kernel
            subplot(2,3,2);     %2. fit of initial parameter guesses
            init_fit_mvg = mvnpdf([x(:) y(:)],mu0,sigma0);
            r_i_mvg = reshape(init_fit_mvg,size(kernel));
            %ifit_norm = rescale(r_i_mvg,[-1 1],[0 255]);
            imagesc(r_i_mvg);
            title('Initial parameters');
            colorbar;
            drawnow;
            subplot(2,3,3);     %3. elliptical 2*sigma error of initial parameters 
            n = length(nonzeros(s_kernel));
            %stdev = 2;
            %pre_factor = 1;
            %[ix,iy] = error_ellipse(mu0,sigma0,stdev,xmean,ymean,n,pre_factor);
            alpha = 0.2;
            i_ellipse = gauss2d_ellipse_byCI(mu0,sigma0,alpha,n);
            %h = plot(xmean+ix,flipud(ymean+iy));
            h = plot(i_ellipse.plot_ellipse(1,:),i_ellipse.plot_ellipse(2,:));
            h.LineWidth = 2.0;
            set(gca,'YDir','reverse');
            xlabel('Space');
            ylabel('Time');
            ylim([1 size(kernel,1)]);
            xlim([1 size(kernel,2)]);
            title('Elliptical fit error of initial parameters')
            drawnow;
            %4. subfield target
            subplot(2,3,4); imagesc(s_kernel); colorbar; title('Sub-field fit target');
            subplot(2,3,5);     %5. live fit
            fit_mvg = mvnpdf([x(:) y(:)],mu,sigma);
            refit_mvg = reshape(fit_mvg,size(kernel));
            %fit_norm = rescale(refit_mvg,scale.*[-1 1],[0 255]);
            imagesc(refit_mvg);
            title(['Live fit. Iteration x: ',num2str(t),' y: ',num2str(u)]);
            colorbar;
            drawnow;
            subplot(2,3,6);     %6. elliptical 2*sigma error of live fit
            %pre_factor_live = 1;
            %[ex,ey] = error_ellipse(mu,sigma,stdev,xmean,ymean,n,pre_factor_live);
            e_ellipse = gauss2d_ellipse_byCI(mu,sigma,alpha,n);
            %hh = plot(xmean+ex,flipud(ymean+ey));
            hh = plot(e_ellipse.plot_ellipse(1,:),e_ellipse.plot_ellipse(2,:));
            hh.LineWidth = 2.0;
            set(gca,'YDir','reverse');
            xlabel('Space');
            ylabel('Time');
            xlim([1 size(kernel,2)]);
            ylim([1 size(kernel,1)]);
            title('Live fit.  Elliptical error');
            drawnow;
            pause(0.1);
        else
        end
        e_array{t,u} = e_ellipse;
    end
end
    
end

