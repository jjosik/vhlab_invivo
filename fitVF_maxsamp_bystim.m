function [ Vth_est_pcoll,Vth_est_ncoll ] = fitVF_maxsamp_bystim( FR_inst,Vm_array,stimvalues,...
    model,h_margins,stimorder,nStim_ON,nStim_OFF,sampling_rate,pref_stim,pref_ind,null_stim,...
    null_ind,tempFrequency,use_flanktest,adapt_bins,anchor_Vth,fit_it,filter_type,save_it,code)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


        
         %%%(OPTIONAL) Select only adjacent directions at peak and null...
         %for V-F analysis (starting with flat criterion of adjacent...
         %directions only; SEE NOTE BELOW)**************
         
         %***||||||||||||||||||||||||||||||||||||||||||||||||||||||
         %***INCORPORATE CODE SELECTING FLANKING DIRECTIONS HERE***
         %***Will require moving up Fourier transform analysis, OR 
         %***moving down VF plot block.
         %***||||||||||||||||||||||||||||||||||||||||||||||||||||||
         
         try
             stim_tf = p.tFrequency;
         catch
             stim_tf = tempFrequency;
         end
            
         criterion = 6;
         
         if use_flanktest == 1,
            %[Vm_ori_th,Vm_ori_accept,sp_ori_th,sp_ori_accept] = ...
            %    run_flankertest(pref_stim,null_stim,stimvalues,...
            %    s_trial_power_coeff,series_power_coeff,freq,s_freq,...
            %    stim_tf,stim_duration,reps_stimmatch,ave_trials,...
            %    criterion);             %*** currently not in effect;
         
         
         else
            all_dirind = 1:length(stimvalues)-1;
            down_ = circshift(stimvalues(1,1:length(stimvalues)-1),1,2);
            up_ = circshift(stimvalues(1,1:length(stimvalues)-1),-1,2);
            down_i = circshift(all_dirind,1,2);
            up_i = circshift(all_dirind,-1,2);
            pref_set = [down_(pref_ind);pref_stim;up_(pref_ind)];
            null_set = [down_(null_ind);null_stim;up_(null_ind)];
            pref_set_index = [down_i(pref_ind);pref_ind;up_i(pref_ind)];
            null_set_index = [down_i(null_ind);null_ind;up_i(null_ind)];
         end
         
         %modify existing Vm traces
         for ii = 1:length(Vm_array),
             new_Vm_array{ii,1} = Vm_array{ii,1}(1:...
                 (length(Vm_array{ii,1})-(round(h_margins*sampling_rate))));
         end
         
         
         array_disp = abs(length(FR_inst{1,1})-length(new_Vm_array{1,1}));
         if length(FR_inst{1,1})==length(new_Vm_array{1,1}),
             coll_pref_Vm = [new_Vm_array{pref_set_index(1,1),1}(1:end-1,1);...
                 new_Vm_array{pref_set_index(2,1),1}(1:end-1,1);...
                 new_Vm_array{pref_set_index(3,1),1}(1:end-1,1)];
             coll_pref_FR = [FR_inst{pref_set_index(1,1),1}(1:end-1,1);...
                 FR_inst{pref_set_index(2,1),1}(1:end-1,1);...
                 FR_inst{pref_set_index(3,1),1}(1:end-1,1)];
             coll_null_Vm = [new_Vm_array{null_set_index(1,1),1}(1:end-1,1);...
                 new_Vm_array{null_set_index(2,1),1}(1:end-1,1);...
                 new_Vm_array{null_set_index(3,1),1}(1:end-1,1)];
             coll_null_FR = [FR_inst{null_set_index(1,1),1}(1:end-1,1);...
                 FR_inst{null_set_index(2,1),1}(1:end-1,1);...
                 FR_inst{null_set_index(3,1),1}(1:end-1,1)];
         else
             for n = 1:length(FR_inst),
                FR_inst{n,1} = reshape(FR_inst{n,1},length(FR_inst{n,1}),1);
             end
             coll_pref_Vm = [new_Vm_array{pref_set_index(1,1),1}(1:end-array_disp,1);...
                 new_Vm_array{pref_set_index(2,1),1}(1:end-array_disp,1);...
                 new_Vm_array{pref_set_index(3,1),1}(1:end-array_disp,1)];
             coll_pref_FR = [FR_inst{pref_set_index(1,1),1};...
                 FR_inst{pref_set_index(2,1),1};...
                 FR_inst{pref_set_index(3,1),1}];
             coll_null_Vm = [new_Vm_array{null_set_index(1,1),1}(1:end-array_disp,1);...
                 new_Vm_array{null_set_index(2,1),1}(1:end-array_disp,1);...
                 new_Vm_array{null_set_index(3,1),1}(1:end-array_disp,1)];
             coll_null_FR = [FR_inst{null_set_index(1,1),1};...
                 FR_inst{null_set_index(2,1),1};...
                 FR_inst{null_set_index(3,1),1}];
         end

         %*****************PREF unbinned*******************************
         stimset = 1;
         z_scale = 1;
         match_PN = 1;  %indicates binning for pref and null plots should be matched
         [ f9,ph,bin_edges ] = VF_densitymap(coll_pref_Vm,coll_null_Vm,...
             coll_pref_FR,adapt_bins,stimset,z_scale,match_PN,...
             pref_set_index,null_set_index);
         cb = colorbar;
         ylabel(cb,'log count of time-step occurances');
         if save_it == 1,
             title(['Maximally-sampled V-F plot'...
                 'for collated PREF direction']);
             saveas(gcf,[pwd filesep strcat(code,'subTH_rawVF_density_PREF.fig')]);
             %close(f9);
         else
         end
         if fit_it == 1,
             %if save_it == 1,
             %    f10 = openfig('subTH_rawVF_density_PREF.fig');
             %    hold on;
             %else
             %end
             sub_elem_pcoll = find(coll_pref_FR==0);
             %Vth_est_pcoll = mode(coll_pref_Vm(sub_elem_pcoll));
             Vth_est_pcoll = mean(coll_pref_Vm(sub_elem_pcoll));
             weight = 1;
             switch model
                 case 0
                     [np_x,pm_out,pcl,pgof,pfitinfo] = vf_powerfit(coll_pref_Vm,...
                         coll_pref_FR,anchor_Vth,Vth_est_pcoll,weight);
                 case 1
                     [np_x,pm_out,pcl,pgof,pfitinfo] = alt_vf_linfit(coll_pref_Vm,...
                         coll_pref_FR,anchor_Vth,Vth_est_pcoll,weight);
             end
             p_z_max = max(max(get(ph,'ZData')));
             np_h = line(np_x,pm_out,p_z_max*ones(length(np_x),1));
             np_h.LineWidth = 3.0;
             np_h.Color = [0 0 0];
             switch model
                 case 0
                     title({'PREF Max-sampled V-F plot ',...
                        'a*rectify(Vm-Vth)^b fit parameters:  Vth = ',...
                        num2str(Vth_est_pcoll),...
                        ' a = ', num2str(pcl.a), ' b = ' num2str(pcl.b)});
                 case 1
                     title({'PREF Max-sampled V-F plot ',...
                        'a*rectify(Vm-Vth) fit parameters: Vth = ',...
                        num2str(Vth_est_pcoll),' a = ', num2str(pcl.a)});
             end
             %mn_line = line([m_pref m_pref],[0 max(pm_out)]);  blocked
             %3/17
             %make local regression smoothing comparison
             %[x,inds] = sort(coll_pref_Vm);
             %yy = smooth(x,coll_pref_FR(inds),5000,'lowess');
             %plot(x,yy,'b-');
             %make gaussian fits across y, smooth means in x bins
             x_bins = [];
             bin_ymean = [];
             bin_ysd = [];
             bin_yse = [];
             bin_count = [];
             for i3 = 2:length(bin_edges{1,1}),
                 if i3 < length(bin_edges{1,1}),
                     bin_inds = find(coll_pref_Vm>=bin_edges{1,1}(1,i3-1)&...
                         coll_pref_Vm<bin_edges{1,1}(1,i3));
                 else
                     bin_inds = find(coll_pref_Vm>=bin_edges{1,1}(1,i3-1)&...
                         coll_pref_Vm<=bin_edges{1,1}(1,i3));
                 end
                 bin_count(end+1,1) = length(bin_inds);
                 if ~isempty(bin_inds),
                     bin_ymean(end+1,1) = nanmean(coll_pref_FR(bin_inds,1));
                     bin_ysd(end+1,1) = std(coll_pref_FR(bin_inds,1));
                     bin_yse(end+1,1) = std(coll_pref_FR(bin_inds,1))./sqrt(length(bin_inds));
                     x_bins = [x_bins;i3-1];
                 else
                 end
             end
             bin_centers = (bin_edges{1,1}(1,2:end)+...
                 bin_edges{1,1}(1,1:end-1))/2;
             scatter3(bin_centers(x_bins),bin_ymean,...
                 (abs(p_z_max+1).*ones(length(bin_ymean),1)));
             for i5 = 1:length(bin_ymean),
                 eb = line([bin_centers(x_bins(i5,1)) bin_centers(x_bins(i5,1))],...
                     [bin_ymean(i5,1)-bin_ysd(i5,1) bin_ymean(i5,1)+bin_ysd(i5,1)]);
                 eb.LineWidth = 2.0;
                 eb.Color = 'r';
             end
             legend('density','fit','1 mV-binned mean Vm',...
                 '1 mV-binned Vm SD','Location','NorthEastOutside');
             hold off;
          else
          end
          if save_it == 1,
              %if filter_type == 2,
              saveas(gcf,[pwd filesep strcat(code,'subTH_maxsampled_VF_density_PREF.fig')]);
                  
              %elseif filter_type == 3,
              %    saveas(gcf,[pwd filesep strcat(code,'boxcar_VF_density_PREF.fig')]);
              %else
              %end
             close(f9);
          else
          end
          %create secondary plot to overlay X-binned pref and null curves
          fel = figure;
          bin_centers = (bin_edges{1,1}(1,2:end)+...
              bin_edges{1,1}(1,1:end-1))/2;
          scatter(bin_centers(x_bins),bin_ymean,'filled','k');
          for i5 = 1:length(bin_ymean),
              eb = line([bin_centers(x_bins(i5,1)) bin_centers(x_bins(i5,1))],...
                  [bin_ymean(i5,1)-bin_yse(i5,1) bin_ymean(i5,1)+bin_yse(i5,1)]);
              eb.LineWidth = 2.0;
              eb.Color = 'r';
          end
          xlabel('Membrane potential (mV)');
          ylabel('Firing rate');
          hold on;
          
          
          
          %plot comparison of FR = 0 and FR > 0 Vm-distributions
          sub_pref_ind = find(coll_pref_FR == 0);
          sub_pref_Vm = coll_pref_Vm(sub_pref_ind);
          prefh = histc(sub_pref_Vm,bin_edges{1,1});
          prefh = prefh(1:end-1);
          fc = figure;
          bar(bin_centers,bin_count);
          grid on;
          grid minor;
          title('Pref response Vm distribution');
          ylabel('Counts');
          xlabel('Potential (mV)');
          hold on;
          bar(bin_centers,prefh,'m');
          legend('All Vm','Subthreshold Vm only',...
              'Location','EastOutside');
          if save_it == 1,
              saveas(gcf,[pwd filesep strcat(code,'PREF_subVm_vs_VmDist','.fig')]);
              close(fc);
          else
          end
          
          %***********************
            %***now NULL unbinned***
            %***********************
          z_scale = 1;
          stimset = null_set_index;
          [ f21,nh,bin_edges_ ] = VF_densitymap(coll_null_Vm,coll_pref_Vm,...
              coll_null_FR,adapt_bins,stimset,z_scale,match_PN,...
              pref_set_index,null_set_index);
          ncb = colorbar;
          ylabel(ncb,'log count of time-steps per joint event bin');
          if save_it == 1,
              title(['Maximally-sampled V-F plot'...
                  'for collated NULL direction']);
              saveas(gcf,[pwd filesep strcat(code,'subTH_rawVF_density_NULL.fig')]);
              %close(f21);
          else
          end
          if fit_it == 1,
              %if save_it == 1,
              %    f22 = openfig('subTH_rawVF_density_NULL.fig');
              %    hold on;
              %else 
              %end
              sub_elem_ncoll = find(coll_null_FR==0);
              %Vth_est_ncoll = mode(coll_null_Vm(sub_elem_ncoll));
              Vth_est_ncoll = mean(coll_null_Vm(sub_elem_ncoll));
              weight = 1;
              switch model
                  case 0
                     [nn_x,nm_out,ncl,ngof,nfitinfo] = vf_powerfit(coll_null_Vm,...
                         coll_null_FR,anchor_Vth,Vth_est_ncoll,weight);
                  case 1
                     [nn_x,nm_out,ncl,ngof,nfitinfo] = alt_vf_linfit(coll_null_Vm,...
                         coll_null_FR,anchor_Vth,Vth_est_ncoll,weight);
              end
              %alt: [nn_x,nm_out,ncl,ngof,nfitinfo] = vf_linfit(coll_null_Vm,coll_null_FR,anchor_Vth,Vth_est_ncoll,weight);
              n_z_max = max(max(get(nh,'ZData')));
              nn_h = line(nn_x,nm_out,n_z_max*ones(length(nn_x),1));
              nn_h.LineWidth = 3.0;
              nn_h.Color = [0 0 0];
              switch model
                  case 0
                      title({'NULL Max-sampled V-F plot ',...
                         'a*rectify(Vm-Vth)^b fit parameters:  Vth ='...
                         num2str(Vth_est_ncoll)...
                         ' a = ' num2str(ncl.a) ' b =' num2str(ncl.b)});
                  case 1
                      title({'NULL Max-sampled V-F plot ',...
                         'a*rectify(Vm-Vth) fit parameters: V_th =',...
                         num2str(Vth_est_ncoll),' a = ', num2str(ncl.a)});
              end
              %mn_line = line([m_null m_null],[0 200]);  ***blocked 3/17
              %make local regression smoothing comparison
              %[xx,inds2] = sort(coll_null_Vm);
              %yyy = smooth(xx,coll_null_FR(inds2),15000,'lowess');
              %plot(xx,yyy,'b-');
              %make gaussian fits across y, smooth means in x bins
              x_bins = [];
              bin_ymean = [];
              bin_ysd = [];
              bin_yse = [];
              bin_count = [];
              for i4 = 2:length(bin_edges{1,1}),
                  if i4 < length(bin_edges{1,1}),
                      bin_inds = find(coll_null_Vm>=bin_edges{1,1}(1,i4-1)&...
                      coll_null_Vm<bin_edges{1,1}(1,i4));
                  else
                      bin_inds = find(coll_null_Vm>=bin_edges{1,1}(1,i4-1)&...
                      coll_null_Vm<=bin_edges{1,1}(1,i4));
                  end
                  bin_count(end+1,1) = length(bin_inds);
                  if ~isempty(bin_inds),
                      bin_ymean(end+1,1) = nanmean(coll_null_FR(bin_inds,1));
                      bin_ysd(end+1,1) = std(coll_null_FR(bin_inds,1));
                      bin_yse(end+1,1) = std(coll_null_FR(bin_inds,1))./sqrt(length(bin_inds));
                      x_bins = [x_bins;i4-1];
                  else
                  end
              end
              bin_centers = (bin_edges{1,1}(1,2:end)+...
              bin_edges{1,1}(1,1:end-1))/2;
              scatter3(bin_centers(x_bins),bin_ymean,...
              (abs(n_z_max+1).*ones(length(bin_ymean),1)));
              for i6 = 1:length(bin_ymean),
                  eb = line([bin_centers(x_bins(i6,1)) bin_centers(x_bins(i6,1))],...
                  [bin_ymean(i6,1)-bin_ysd(i6,1) bin_ymean(i6,1)+bin_ysd(i6,1)]);
                  eb.LineWidth = 2.0;
                  eb.Color = 'b';
              end
            legend('density','fit','1 mV-binned mean Vm',...
                '1 mV-binned Vm SD','Location','NorthEastOutside');
            hold off;
        else
        end
        if save_it == 1,
            %if filter_type == 2,
            saveas(gcf,[pwd filesep strcat(code,'_subTH_maxsampled_VF_density_NULL.fig')]);
                
            %elseif filter_type == 3,
            %      saveas(gcf,[pwd filesep strcat(code,'boxcar_VF_density_NULL.fig')]);
            %else
            %end
            close(f21);
        else
        end
        %complete secondary plot
        figure(fel);
        bin_centers = (bin_edges{1,1}(1,2:end)+...
            bin_edges{1,1}(1,1:end-1))/2;
        scatter(bin_centers(x_bins),bin_ymean,'filled','k');
        for i6 = 1:length(bin_ymean),
            eb = line([bin_centers(x_bins(i6,1)) bin_centers(x_bins(i6,1))],...
                [bin_ymean(i6,1)-bin_yse(i6,1) bin_ymean(i6,1)+bin_yse(i6,1)]);
            eb.LineWidth = 2.0;
            eb.Color = 'b';
        end
        title('X-binned Firing rate means and errors, PREF and NULL');
        hold off;
        if save_it == 1,
            saveas(gcf,[pwd filesep strcat(code,'Xbin_VF_overlay.fig')]);
            close(fel);
        else
        end
        
        
         %plot comparison of FR = 0 and FR > 0 Vm-distributions
        sub_null_ind = find(coll_null_FR == 0);
        sub_null_Vm = coll_null_Vm(sub_null_ind);
        nullh = histc(sub_null_Vm,bin_edges{1,1});
        nullh = nullh(1:end-1);
        ffn = figure;
        bar(bin_centers,bin_count);
        grid on;
        grid minor;
        title('Null response Vm distribution');
        ylabel('Counts');
        xlabel('Potential (mV)');
        hold on;
        bar(bin_centers,nullh,'c');
        legend('All Vm','Subthreshold Vm only','Location','NorthEast');
        if save_it == 1,
            saveas(gcf,[pwd filesep 'PREF_subVm_vs_VmDist','.fig']);
            close(ffn);
        else
        end
        
        if save_it == 1,
            vffilename = strcat(code,'_fitdata.mat');
            save(vffilename);
        else
        end



end

