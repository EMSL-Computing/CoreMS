
from numpy.fft import rfft
from pywt import dwt 
from scipy.spatial.distance import cosine, jaccard, euclidean, cityblock
from scipy.stats import pearsonr, spearmanr, kendalltau
from sklearn.metrics.pairwise import cosine_similarity
from numpy import power, dot, absolute, subtract, intersect1d, where, average, corrcoef
from numpy.linalg import norm
from pandas import DataFrame

def weighted_cosine_correlation(  ms_mz_abun_dict, ref_obj, a=0.5, b=1.3):
    
    # create dict['mz'] = abundance, for experimental data
    #ms_mz_abun_dict = mass_spec.mz_abun_dict

    # weight exp data

    exp_abun = list(ms_mz_abun_dict.values())
    exp_mz = list(ms_mz_abun_dict.keys())

    xc = power(exp_abun, a) *  power(exp_mz, b) 
    
    # track back to individual mz
    weighted_exp_dict = dict(zip(ms_mz_abun_dict.keys(), xc))

    # weight ref data
    yc = power(ref_obj.get("abundance"), a) *  power(ref_obj.get("mz"), b) 
    
    ref_mz_abun_dict = dict(zip(ref_obj.get("mz"), yc))

    #parse to dataframe, easier to zerofill and tranpose
    df = DataFrame([weighted_exp_dict, ref_mz_abun_dict])

    # fill missing mz with weight {abun**a}{m/z**b} to 0
    df.fillna(0, inplace=True)
    
    #calculate cosine correlation, 
    x = df.T[0].values
    y = df.T[1].values

    #correlation = (1 - cosine(x, y))
    
    correlation = dot(x, y)/(norm(x)*norm(y))

    return correlation

def cosine_correlation(  ms_mz_abun_dict, ref_obj):

    # create dict['mz'] = abundance, for experimental data
    #ms_mz_abun_dict = mass_spec.mz_abun_dict

    # create dict['mz'] = abundance, for experimental data
    ref_mz_abun_dict = dict(zip(ref_obj.get("mz"), ref_obj.get("abundance")))

    #print(ref_mz_abun_dict)
    #print(ms_mz_abun_dict)

    #parse to dataframe, easier to zerofill and tranpose
    df = DataFrame([ms_mz_abun_dict, ref_mz_abun_dict])

    # fill missing mz with abundance 0
    df.fillna(0, inplace=True)
    
    #calculate cosine correlation, 
    x = df.T[0].values
    y = df.T[1].values
    
    #correlation = (1 - cosine(x, y))
    
    correlation = dot(x, y)/(norm(x)*norm(y))

    return correlation

def stein_scott( ms_mz_abun_dict, ref_obj):
 
    #ms_mz_abun_dict = mass_spec.mz_abun_dict
 
    exp_abun = list(ms_mz_abun_dict.values())
    exp_mz = list(ms_mz_abun_dict.keys())
 
    ref_abun = ref_obj.get("abundance")
    ref_mz = ref_obj.get("mz")
 
    # important: I assume ref_mz and ref_abun are in order, and one-to-one; this needs to be be verified
    ref_mz_abun_dict = dict(zip(ref_mz, ref_abun))
 
    # filter out the mass values that have zero intensities in exp_mz
    exp_mz_filtered = set([k for k in exp_mz if ms_mz_abun_dict[k] != 0])
 
    # filter out the mass values that have zero intensities in ref_mz
    ref_mz_filtered = set([k for k in ref_mz if ref_mz_abun_dict[k] != 0])
 
    # find the intersection/common mass values of both ref and exp, and sort them
    common_mz_values = sorted(list(exp_mz_filtered.intersection(ref_mz_filtered)))
    
    # find the number of common mass values (after filtering 0s)
    n_x_y = len(common_mz_values)
    
    if n_x_y == 0: return 0
    
    # count number of non-zero abundance/peak intensity values
    n_x = sum(a != 0 for a in exp_abun)
 
    s_r_x_y = 0
    
    a, b = 1, 0
    for i in range(1,n_x_y):
        
        current_value = common_mz_values[i]
        previous_value = common_mz_values[i-1]
 
        y_i = ref_mz_abun_dict[current_value]
        y_i_minus1 = ref_mz_abun_dict[previous_value]
        
        lc_current = power(y_i, a) *  power(current_value, b)
        lc_previous = power(y_i_minus1, a) *  power(previous_value, b)
        
        x_i = ms_mz_abun_dict[current_value]
        x_i_minus1 = ms_mz_abun_dict[previous_value]
        
        uc_current = power(x_i, a) *  power(current_value, b)
        uc_previous = power(x_i_minus1, a) *  power(previous_value, b)

        T1 = lc_current/lc_previous
        
        T2 = uc_previous/uc_current

        temp_computation = T1 * T2
        
        n = 0
        if temp_computation <= 1:
            n = 1
        else:
            n = -1
        s_r_x_y = s_r_x_y + power(temp_computation,n)
 
    # finish the calculation of S_R(X,Y)
    
    s_r_x_y = s_r_x_y / n_x_y
    # using the existing weighted_cosine_correlation function to get S_WC(X,Y)
    s_wc_x_y = weighted_cosine_correlation(ms_mz_abun_dict, ref_obj, a=0.5, b =3 )
 
    # final step
    s_ss_x_y = ( (n_x * s_wc_x_y) + (n_x_y * s_r_x_y) )/ (n_x + n_x_y)
 
    return s_ss_x_y
 

def pearson_correlation(  ms_mz_abun_dict, ref_obj):

    # create dict['mz'] = abundance, for experimental data
    #ms_mz_abun_dict = mass_spec.mz_abun_dict

    # create dict['mz'] = abundance, for experimental data
    ref_mz_abun_dict = dict(zip(ref_obj.get("mz"), ref_obj.get("abundance")))

    # parse to dataframe, easier to zerofill and tranpose
    df = DataFrame([ms_mz_abun_dict, ref_mz_abun_dict])

    # fill missing mz with abundance 0
    df.fillna(0, inplace=True)
    
    # calculate Pearson correlation
    correlation = pearsonr(df.T[0], df.T[1])

    return correlation[0]

def spearman_correlation(  ms_mz_abun_dict, ref_obj):

    # create dict['mz'] = abundance, for experimental data
    #ms_mz_abun_dict = mass_spec.mz_abun_dict

    # create dict['mz'] = abundance, for experimental data
    ref_mz_abun_dict = dict(zip(ref_obj.get("mz"), ref_obj.get("abundance")))

    # parse to dataframe, easier to zerofill and tranpose
    df = DataFrame([ms_mz_abun_dict, ref_mz_abun_dict])

    # fill missing mz with abundance 0
    df.fillna(0, inplace=True)
    
    # calculate Spearman correlation
    ### TODO - Check axis
    correlation = spearmanr(df.T[0], df.T[1], axis=0)

    return correlation[0]

def kendall_tau(  ms_mz_abun_dict, ref_obj):

    # create dict['mz'] = abundance, for experimental data
    #ms_mz_abun_dict = mass_spec.mz_abun_dict

    # create dict['mz'] = abundance, for experimental data
    ref_mz_abun_dict = dict(zip(ref_obj.get("mz"), ref_obj.get("abundance")))

    # parse to dataframe, easier to zerofill and tranpose
    df = DataFrame([ms_mz_abun_dict, ref_mz_abun_dict])

    # fill missing mz with abundance 0
    df.fillna(0, inplace=True)
    
    # calculate Kendall's tau
    correlation = kendalltau(df.T[0], df.T[1])

    return correlation[0]

def euclidean_distance(  ms_mz_abun_dict, ref_obj):

    def euclidean_distance_manual(qlist,rlist):

        T1=sum(subtract(qlist,rlist)**2)

        T2=sum((qlist)**2)

        return (1+T1/T2)**(-1)
        
    # create dict['mz'] = abundance, for experimental data
    #ms_mz_abun_dict = mass_spec.mz_abun_dict

    # create dict['mz'] = abundance, for experimental data
    ref_mz_abun_dict = dict(zip(ref_obj.get("mz"), ref_obj.get("abundance")))

    # parse to dataframe, easier to zerofill and tranpose
    df = DataFrame([ms_mz_abun_dict, ref_mz_abun_dict])

    # fill missing mz with abundance 0
    df.fillna(0, inplace=True)
    
    # calculate Pearson correlation
    #correlation = euclidean(df.T[0], df.T[1])

    correlation = euclidean_distance_manual(df.T[0], df.T[1])
    return correlation

def manhattan_distance(  ms_mz_abun_dict, ref_obj):
    
    def mann_distance_manual(qlist,rlist):
    
        T1=sum(absolute(subtract(qlist,rlist)))
        T2=sum(qlist)
        return (1+T1/T2)**(-1)

    # create dict['mz'] = abundance, for experimental data
    #ms_mz_abun_dict = mass_spec.mz_abun_dict

    # create dict['mz'] = abundance, for experimental data
    ref_mz_abun_dict = dict(zip(ref_obj.get("mz"), ref_obj.get("abundance")))

    # parse to dataframe, easier to zerofill and tranpose
    df = DataFrame([ms_mz_abun_dict, ref_mz_abun_dict])

    # fill missing mz with abundance 0
    df.fillna(0, inplace=True)
    
    # calculate manhattan correlation
    #correlation = cityblock(df.T[0], df.T[1])
    correlation = mann_distance_manual(df.T[0], df.T[1])
    
    return correlation

def dft_correlation(ms_mz_abun_dict, ref_obj):

    #ms_mz_abun_dict = mass_spec.mz_abun_dict

    exp_abun = list(ms_mz_abun_dict.values())
    exp_mz = list(ms_mz_abun_dict.keys())

    ref_abun = ref_obj.get("abundance")
    ref_mz = ref_obj.get("mz")

    # important: I assume ref_mz and ref_abun are in order, and one-to-one; this needs to be be verified
    ref_mz_abun_dict = dict(zip(ref_mz, ref_abun))

    # filter out the mass values that have zero intensities in exp_mz
    exp_mz_filtered = set([k for k in exp_mz if ms_mz_abun_dict[k] != 0])

    # filter out the mass values that have zero intensities in ref_mz
    ref_mz_filtered = set([k for k in ref_mz if ref_mz_abun_dict[k] != 0])

    # find the intersection/common mass values of both ref and exp, and sort them
    
    common_mz_values = exp_mz_filtered.intersection(ref_mz_filtered)
    
    # find the number of common mass values (after filtering 0s)
    n_x_y = len(common_mz_values)
    
    if n_x_y == 0: return 0

    # count number of non-zero abundance/peak intensity values
    n_x = sum(a != 0 for a in exp_abun)

    #parse to dataframe, easier to zerofill and tranpose
    df = DataFrame([ms_mz_abun_dict, ref_mz_abun_dict])

    # fill missing mz with abundance 0
    df.fillna(0, inplace=True)
    
    #calculate cosine correlation, 
    x = df.T[0].values
    y = df.T[1].values

    # get the Fourier transform of x and y
    x_dft = rfft(x).real
    y_dft = rfft(y).real

    s_dft_xy = dot(x_dft, y_dft)/(norm(x_dft)*norm(y_dft))

    # using the existing weighted_cosine_correlation function to get S_WC(X,Y)
    s_wc_x_y = weighted_cosine_correlation(ms_mz_abun_dict, ref_obj)

    # final step
    s_dft = (n_x * s_wc_x_y + n_x_y * s_dft_xy) / (n_x + n_x_y)

    return s_dft

def jaccard_distance(  ms_mz_abun_dict, ref_obj):
    
    def jaccard_similarity(list1, list2):
        intersection = len(list(set(list1).intersection(list2)))
        union = (len(list1) + len(list2)) - intersection
        return float(intersection) / union
    
    # create dict['mz'] = abundance, for experimental data
    #ms_mz_abun_dict = mass_spec.mz_abun_dict

    # create dict['mz'] = abundance, for experimental data
    ref_mz_abun_dict = dict(zip(ref_obj.get("mz"), ref_obj.get("abundance")))

    # parse to dataframe, easier to zerofilland tranpose
    df = DataFrame([ms_mz_abun_dict, ref_mz_abun_dict])

    # fill missing mz with abundance 0
    df.fillna(0, inplace=True)
    
    # calculate jaccard correlation
    #correlation = jaccard(df.T[0], df.T[1])
    correlation = jaccard_similarity(df.T[0], df.T[1])
    return correlation

def dwt_correlation(  ms_mz_abun_dict, ref_obj):
    #Need to have pywavelets installed
    #ms_mz_abun_dict = mass_spec.mz_abun_dict

    exp_abun = list(ms_mz_abun_dict.values())
    exp_mz = list(ms_mz_abun_dict.keys())

    ref_abun = ref_obj.get("abundance")
    ref_mz = ref_obj.get("mz")

    # important: I assume ref_mz and ref_abun are in order, and one-to-one; this needs to be be verified
    ref_mz_abun_dict = dict(zip(ref_mz, ref_abun))

    # filter out the mass values that have zero intensities in exp_mz
    exp_mz_filtered = set([k for k in exp_mz if ms_mz_abun_dict[k] != 0])

    # filter out the mass values that have zero intensities in ref_mz
    ref_mz_filtered = set([k for k in ref_mz if ref_mz_abun_dict[k] != 0])

    # find the intersection/common mass values of both ref and exp, and sort them
    common_mz_values = sorted(list(exp_mz_filtered.intersection(ref_mz_filtered)))
    
    # find the number of common mass values (after filtering 0s)
    n_x_y = len(common_mz_values)

    if n_x_y == 0: return 0

    # count number of non-zero abundance/peak intensity values
    n_x = sum(a != 0 for a in exp_abun)

    #parse to dataframe, easier to zerofill and tranpose
    df = DataFrame([ms_mz_abun_dict, ref_mz_abun_dict])

    # fill missing mz with abundance 0
    df.fillna(0, inplace=True)
    
    #calculate cosine correlation, 
    x = df.T[0].values
    y = df.T[1].values

    #Make x and y into an array
    x_a=list(x)
    y_a=list(y)

    # get the wavelet transform of x and y (Daubechies with a filter length of 4. Asymmetric. pywavelets function)
    #Will only use the detail dwt (dwtDd
    (x_dwtA,x_dwtD) = dwt(x_a,'db2')
    (y_dwtA,y_dwtD) = dwt(y_a,'db2')

    s_dwt_xy = dot(x_dwtD, y_dwtD)/(norm(x_dwtD)*norm(y_dwtD))

    # using the existing weighted_cosine_correlation function to get S_WC(X,Y)
    s_wc_x_y = weighted_cosine_correlation(ms_mz_abun_dict, ref_obj)

    # final step
    s_dwt = (n_x * s_wc_x_y + n_x_y * s_dwt_xy) / (n_x + n_x_y)

    return s_dwt