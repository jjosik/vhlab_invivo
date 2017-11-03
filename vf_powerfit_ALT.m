function [ n_x,m_out,cl,gof,fitinfo ] = vf_powerfit_ALT( Vm,FR_inst,anchor_Vth,Vth_est,weight )

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
                'Lower',[0,Vth_est,0.1],...
                'Upper',[500,Vth_est,8],...
                'Startpoint',[5,Vth_est,1]);
        else
            s = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',[0,Vth_est,0.1],...
                'Upper',[500,Vth_est,8],...
                'Startpoint',[5,Vth_est,1]);
        end
        ft = fittype('a.*rectify(x-Vth_est).^b','options',s);
        [cl,gof,fitinfo] = fit(Vm,FR_inst,ft);
        n_x = [floor(min(Vm))+0.5:1:floor(max(Vm))+0.5];
        m_out = cl.a.*(n_x-Vth_est).^cl.b;
        m_out = heaviside(n_x-Vth_est).*m_out;
    else
        if weight == 1,
            s = fitoptions('Method','NonlinearLeastSquares',...
                'Weight',repmat(weight,length(FR_inst),1),...
                'Lower',[0,min(Vm),0.1],...
                'Upper',[500,max(Vm),8],...
                'Startpoint',[5,min(Vm),1]);
        else
            s = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',[0,min(Vm),0.1],...
                'Upper',[500,max(Vm),8],...
                'Startpoint',[5,min(Vm),1]);
        end
        ft = fittype('a.*rectify(x-b).^c','options',s);
        [cl,gof,fitinfo] = fit(Vm,FR_inst,ft);
        n_x = [floor(min(Vm))+0.5:1:floor(max(Vm))+0.5];
        m_out = cl.a.*(n_x-cl.b).^cl.c;
        m_out = heaviside(n_x-cl.b).*m_out;
    end



end



