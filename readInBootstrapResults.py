import numpy as np
import matplotlib.pyplot as plt
from cantrips import union
import manageMCMCResults as mMCMCR
from AstronomicalParameterArchive import AstronomicalParameterArchive

def unpackOneResultsSet(results):
    return 0

def plotBootstrapResults(bootstrap_results, trim_portions = [0.05, 0.95], max_bin_display_scaling = 1.01, n_bins = 20, save = 0, show = 0, save_dir = '/Users/sasha/Documents/Harvard/physics/randall/randomizationResults/', save_file_name = 'bootstrap_results.png'):
    plt.rc('font', family='serif')
    plt.rc('text', usetex=True)
    #param_display_vals = mMCMCR.getParamDisplayVals()
    astro_arch = AstronomicalParameterArchive()
    deg_to_rad = astro_arch.getDegToRad()
    param_display_labels = mMCMCR.getDisplayLabels()
    halo_types = list(bootstrap_results.keys())
    all_params = union([list(bootstrap_results[halo_type].keys()) for halo_type in halo_types])
    params_to_use = ['M', 'rs', 'el', 'phi', 'h_x_center', 'h_z_center']
    halo_titles = {'nfw_obl':'Oblate NFW DM Halo', 'nfw_pro':'Prolate NFW DM Halo',
                   'cored_obl':'Oblate Acored DM Halo', 'cored_pro':'Prolate Acored DM Halo',
                   'burkert_obl':'Oblate Burkert DM Halo', 'burkert_pro':'Prolate Burkert DM Halo'}
    param_scalings = [1.0, 1.0, 1.0, 1.0 / deg_to_rad, 1.0, 1.0]
    n_params = len(params_to_use)
    n_halo_types = len(halo_types)
    f, big_axarr = plt.subplots(figsize = (16, 16), nrows = len(halo_types), ncols = 1, sharey = True)
    f.subplots_adjust(wspace = 0.05, hspace = 0.5)

    for i in range(len(halo_types)):
        halo_type = halo_types[i]
        title = (halo_titles[halo_type] if halo_type in halo_titles.keys() else halo_type)
        print ('title = ' + str(title))
        big_axarr[i].set_title(title)
        big_axarr[i].tick_params(labelcolor=(1.,1.,1., 0.0), top='off', bottom='off', left='off', right='off')
        big_axarr[i]._frameon = False

    for i in range(len(halo_types)):
        halo_type = halo_types[i]

        print ('Computing results for halo type ' + str(halo_type))
        single_bootstrap_results = bootstrap_results[halo_type]
        max_bin_height = 0.0
        new_axes = [f.add_subplot(len(halo_types), len(params_to_use) , i * len(params_to_use) + j + 1) for j in range(len(params_to_use))]
        for j in range(len(params_to_use)):
            ax = new_axes[j]
            param = params_to_use[j]
            param_scaling = param_scalings[j]
            if param in list(single_bootstrap_results.keys()):
                single_param_results = np.array(single_bootstrap_results[param]) * param_scaling
                single_param_results = sorted(single_param_results)
                n_results = len(single_param_results)
                single_hist = ax.hist(single_param_results, bins = n_bins, histtype = 'step', color = 'black')
                max_bin_height = max([max_bin_height] + list(single_hist[0]))
                #if i == n_halo_types - 1:
                if param in list(param_display_labels.keys()):
                    ax.set_xlabel(param_display_labels[param])
                else:
                    ax.set_xlabel(param)
            else:
                print ('param ' + str(param) + ' not in results for halo_type ' + str(halo_type) )
            if j > 0:
                ax.set_yticklabels([])
            else:
                ax.set_ylabel('Number of Bootstraps')
        for ax in new_axes: ax.set_ylim([0.0, max_bin_height * max_bin_display_scaling])
    if save:
        plt.savefig(save_dir + save_file_name)
    if show:
        plt.show()
    plt.close('all')
    return 0

def measureBootstrapRatios(bootstrap_results, likelihood_key = 'likelihood'):
    likelihood_sets = {}
    halo_types = list(bootstrap_results.keys())
    for halo_type in halo_types:
        likelihood_sets[halo_type] = bootstrap_results[halo_type]['likelihood'][:]
    halo_ratios = {}
    for i in range(len(halo_types)):
        halo_type1 = halo_types[i]
        for j in range(len(halo_types)):
            halo_type2 = halo_types[j]
            if i != j:
                print ('[halo_type1, halo_type2] = ' + str([halo_type1, halo_type2]))
                likelihood_ratios = np.exp(np.array(bootstrap_results[halo_type1][likelihood_key]) - np.array(bootstrap_results[halo_type2][likelihood_key]))
                log10_likelihood_ratios = np.log10(likelihood_ratios)
                n_larger_rats = len([rat for rat in likelihood_ratios if rat > 1.0])
                median_log10_likelihood_rat = np.median(log10_likelihood_ratios)
                mean_log10_likelihood_rat = np.mean(log10_likelihood_ratios)
                std_log10_likelihood_rat = np.std(log10_likelihood_ratios)
                halo_ratios[halo_type1 + '/' + halo_type2] = [likelihood_ratios, n_larger_rats, median_log10_likelihood_rat, mean_log10_likelihood_rat, std_log10_likelihood_rat]

    return halo_ratios



#bootstrap results are of form { 'halo_type1': {'param1': [res1, res2...], 'param2': [res1, res2...]},
#                                'halo_type2': {'param1': [res1, res2...], 'param2': [res1, res2...]}, ... }
def measureBootstrapResults(bootstrap_results, trim_portions = [0.05, 0.95]):

    #astro_arch = AstronomicalParameterArchive()
    #deg_to_rad = astro_arch.getDegToRad()
    for halo_type in bootstrap_results.keys():
        print ('Computing results for halo type ' + str(halo_type))
        single_bootstrap_results = bootstrap_results[halo_type]
        for param in single_bootstrap_results.keys():
            single_param_results = single_bootstrap_results[param]
            single_param_results = sorted(single_param_results)
            n_results = len(single_param_results)
            single_param_results = single_param_results[int(trim_portions[0] * n_results):int(trim_portions[1] * n_results)]
            mean = np.mean(single_param_results)
            median = np.median(single_param_results)
            std = np.std(single_param_results)
            print ('For [halo_type, param] ' + str([halo_type, param]) + ': [median, mean, std] = ' + str([median, mean, std]))

    return 0


def readInBootstrapResults(bootstrap_files, params_to_read_in = None, elliptical = 0): #params_to_read_in = ['el', 'phi', 'h_x_center', 'rs', 'M', 'h_z_center']):
    if params_to_read_in == None:
        if elliptical:
            params_to_read_in = ['M', 'rs', 'el','phi', 'sigsqr_RR_MP', 'sigsqr_RR_IM', 'sigsqr_RR_MR', 'omega_phi_MP', 'omega_phi_IM', 'omega_phi_MR']
        else:
            params_to_read_in = ['M', 'rs', 'h_x_center','h_z_center', 'sigsqr_rr_0_MP', 'sigsqr_rr_0_IM', 'sigsqr_rr_0_MR', ]

    n_items_per_halo = 3
    name_index = 0
    fitted_params_index = 1
    likelihood_index = 2
    print ('Hi')
    if elliptical:
        first_results_set = np.load(bootstrap_files[0], encoding = 'latin1', allow_pickle = True)[0]
    else:
        first_results_set = np.load(bootstrap_files[0], encoding = 'latin1', allow_pickle = True)
    print('first_results_set =  ' + str(first_results_set))
    #first_results_set = np.load(bootstrap_files[0], allow_pickle = True)
    #first_results = first_results_set[0]

    n_halos = len(first_results_set) // n_items_per_halo
    print ('n_halos = ' + str(n_halos))
    halo_types = [first_results_set[name_index + n_items_per_halo * i] for i in range(n_halos) ]
    print ('halo_types = '  + str(halo_types))

    param_dict = {}
    for param in params_to_read_in:
        param_dict[param] = []
    param_dict['likelihood'] = []

    dicts_by_halo = {}
    for halo in halo_types:
        dicts_by_halo[halo] = {}
        for param in params_to_read_in:
            dicts_by_halo[halo][param] = []
        dicts_by_halo[halo]['likelihood'] = []

    if elliptical:
        all_results = [np.load(load_file, encoding = 'latin1', allow_pickle = True)[0] for load_file in bootstrap_files]
    else:
        all_results = [np.load(load_file, encoding = 'latin1', allow_pickle = True) for load_file in bootstrap_files]

    print ('all_results = ' + str(all_results))
    for run_results in all_results:
        print ('type(run_results) = ' + str(type(run_results) ))
        print ('run_results = ' + str(run_results))
        run_results = {run_results[i * 3 + name_index]:[run_results[i * 3 + fitted_params_index], run_results[i * 3 + likelihood_index]] for i in range(len(run_results) // 3)}
        print ('run_results = ' + str(run_results))
        print('run_results.keys() = ' + str(run_results.keys()))
        #for single_results_set in run_results:
        #    print ('single_results_set = ' + str(single_results_set))
        for halo_key in list(run_results.keys()):
                #halo = single_results_set[name_index + i * n_items_per_halo]
                #fitted_params = single_results_set[fitted_params_index + i * n_items_per_halo]
                #likelihood = single_results_set[likelihood_index + i * n_items_per_halo]
                fitted_params = run_results[halo_key][0]
                likelihood = run_results[halo_key][1]
                for param in params_to_read_in:
                    print ('[param, halo_key, likelihood] = ' + str([param, halo_key, likelihood] ))
                    dicts_by_halo[halo_key][param] = dicts_by_halo[halo_key][param] + [fitted_params[param] ]
                dicts_by_halo[halo_key]['likelihood'] = dicts_by_halo[halo_key]['likelihood'] + [likelihood]
                print ('dicts_by_halo["halo_key"]["likelihood"] = ' + str(dicts_by_halo[halo_key]["likelihood"]))

    return dicts_by_halo
