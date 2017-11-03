function [ n_x,m_out,cl,gof,fitinfo ] = alt_vf_linfit( Vm,FR_inst,anchor_Vth,Vth_est,weight )

% only difference between this version and "vf_linfit.m" is the use of
% "linearized" log-exp fit function

%NOTE - observe occasional error where the order of the parameters get
%shuffled; the current ordering in fitoptions reflects the presumption that
%the ordering is implicitly determined by the order of operations in the
%model equation.  Seems to work for now but beware of potential for future errors.
FR_inst = reshape(FR_inst,length(FR_inst),1);
old = sympref('HeavisideAtOrigin', 0);
    if anchor_Vth == 1,
        if weight == 1,
            s = fitoptions('Method','NonlinearLeastSquares',...
                'Weight',repmat(weight,length(FR_inst),1),...
                'Lower',[Vth_est,0],...
                'Upper',[Vth_est,500],...
                'Startpoint',[Vth_est, 5]);
        else
            s = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',[Vth_est,0],...
                'Upper',[Vth_est,500],...
                'Startpoint',[Vth_est, 5]);
           
        end
        ft = fittype('a.*rectify(log(1+exp(x-Vth_est)))','options',s);
        [cl,gof,fitinfo] = fit(Vm,FR_inst,ft);
        n_x = [min(Vm):1:max(Vm)];
        m_out = cl.a.*(log(1+exp(n_x-Vth_est)));
        m_out = heaviside(log(1+exp(n_x-Vth_est))).*m_out;
    else
        if weight == 1,
            s = fitoptions('Method','NonlinearLeastSquares',...
                'Weight',repmat(weight,length(FR_inst),1),...
                'Lower',[min(Vm),0],...
                'Upper',[max(Vm),500],...
                'Startpoint',[min(Vm),5]);
            
        else
            s = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',[min(Vm),0],...
                'Upper',[max(Vm),500],...
                'Startpoint',[min(Vm),5]);
           
        end
        ft = fittype('a.*rectify(log(1+exp(x-b))','options',s);
        [cl,gof,fitinfo] = fit(Vm,FR_inst,ft);
        n_x = [min(Vm):1:max(Vm)];
        m_out = cl.a.*(log(1+exp(n_x-cl.b)));
        m_out = heaviside(log(1+exp(n_x-Vth_est))).*m_out;
    end



end



