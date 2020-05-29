
from numpy.fft import rfft
from pywt import dwt 
from scipy.spatial.distance import cosine, jaccard, euclidean, cityblock
from scipy.stats import pearsonr, spearmanr, kendalltau
from sklearn.metrics.pairwise import cosine_similarity
from numpy import power, dot, absolute, subtract, intersect1d, where, average, corrcoef
from numpy.linalg import norm
from pandas import DataFrame

class SpectralSimilarity():

    def __init__(self, ms_mz_abun_dict, ref_obj):
        
        self.ms_mz_abun_dict = ms_mz_abun_dict
        self.ref_obj = ref_obj
        
        self.exp_abun = list(self.ms_mz_abun_dict.values())
        self.exp_mz = list(self.ms_mz_abun_dict.keys())

        self.ref_mz = self.ref_obj.get("mz")
        self.ref_abun = self.ref_obj.get("abundance")

        self.ref_mz_abun_dict = dict(zip(self.ref_mz , self.ref_abun))

        #parse to dataframe, easier to zerofill and tranpose
        df = DataFrame([self.ms_mz_abun_dict, self.ref_mz_abun_dict])

        # fill missing mz with abundance 0
        df.fillna(0, inplace=True)
        
        #calculate cosine correlation, 
        x = df.T[0].values
        y = df.T[1].values
        
        self.zero_filled_u_l = (x,y)

        # filter out the mass values that have zero intensities in self.exp_abun
        exp_mz_filtered = set([k for k in self.exp_mz if self.ms_mz_abun_dict[k] != 0])
    
        # filter out the mass values that have zero intensities in self.ref_mz
        self.ref_mz_filtered = set([k for k in self.ref_mz if self.ref_mz_abun_dict[k] != 0])
    
        # find the intersection/common mass values of both ref and exp, and sort them
        self.common_mz_values = sorted(list(exp_mz_filtered.intersection(self.ref_mz_filtered)))
        
        # find the number of common mass values (after filtering 0s)
        self.n_x_y = len(self.common_mz_values)


    def weighted_cosine_correlation(self, a=0.5, b=1.3):
        
        # create dict['mz'] = abundance, for experimental data
        #ms_mz_abun_dict = mass_spec.mz_abun_dict

        # weight exp data

        xc = power(self.exp_abun, a) *  power(self.exp_abun, b) 
        
        # track back to individual mz
        weighted_exp_dict = dict(zip(self.ms_mz_abun_dict.keys(), xc))

        # weight ref data
        yc = power(self.ref_obj.get("abundance"), a) *  power(self.ref_obj.get("mz"), b) 
        
        ref_mz_abun_dict = dict(zip(self.ref_obj.get("mz"), yc))

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

    def cosine_correlation(self):

        #calculate cosine correlation, 
        x = self.zero_filled_u_l[0]
        y = self.zero_filled_u_l[1]
        
        #correlation = (1 - cosine(x, y))
        
        correlation = dot(x, y)/(norm(x)*norm(y))

        return correlation

    def stein_scott(self):
    
        if self.n_x_y == 0: return 0
        
        # count number of non-zero abundance/peak intensity values
        n_x = sum(a != 0 for a in self.exp_abun)
    
        s_r_x_y = 0
        
        a, b = 1, 0

        for i in range(1,self.n_x_y):
            
            current_value = self.common_mz_values[i]
            previous_value = self.common_mz_values[i-1]
    
            y_i = self.ref_mz_abun_dict[current_value]
            y_i_minus1 = self.ref_mz_abun_dict[previous_value]
            
            lc_current = power(y_i, a) *  power(current_value, b)
            lc_previous = power(y_i_minus1, a) *  power(previous_value, b)
            
            x_i = self.ms_mz_abun_dict[current_value]
            x_i_minus1 = self.ms_mz_abun_dict[previous_value]
            
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
        
        s_r_x_y = s_r_x_y / self.n_x_y
        # using the existing weighted_cosine_correlation function to get S_WC(X,Y)
        s_wc_x_y = self.weighted_cosine_correlation(a=0.5, b =3 )
    
        # final step
        s_ss_x_y = ( (n_x * s_wc_x_y) + (self.n_x_y * s_r_x_y) )/ (n_x + self.n_x_y)
    
        return s_ss_x_y
    

    def pearson_correlation(self,):

        correlation = pearsonr(self.zero_filled_u_l[0], self.zero_filled_u_l[1])

        return correlation[0]

    def spearman_correlation(self):

        # calculate Spearman correlation
        ### TODO - Check axis
        correlation = spearmanr(self.zero_filled_u_l[0], self.zero_filled_u_l[1], axis=0)

        return correlation[0]

    def kendall_tau(self):

        # create dict['mz'] = abundance, for experimental data
        #self.ms_mz_abun_dict = mass_spec.mz_abun_dict

        # create dict['mz'] = abundance, for experimental data
        
        # calculate Kendall's tau
        correlation = kendalltau(self.zero_filled_u_l[0], self.zero_filled_u_l[1])

        return correlation[0]

    def euclidean_distance(self):

        def euclidean_distance_manual(qlist,rlist):

            T1=sum(subtract(qlist,rlist)**2)

            T2=sum((qlist)**2)

            return (1+T1/T2)**(-1)
        
        correlation = euclidean_distance_manual(self.zero_filled_u_l[0], self.zero_filled_u_l[1])
        return correlation

    def manhattan_distance(self):
        
        def mann_distance_manual(qlist,rlist):
        
            T1=sum(absolute(subtract(qlist,rlist)))
            T2=sum(qlist)
            return (1+T1/T2)**(-1)

        # calculate manhattan correlation
        #correlation = cityblock(df.T[0], df.T[1])
        correlation = mann_distance_manual(self.zero_filled_u_l[0], self.zero_filled_u_l[1])
        
        return correlation

    def dft_correlation(self):

        if self.n_x_y == 0: return 0

        # count number of non-zero abundance/peak intensity values
        n_x = sum(a != 0 for a in self.exp_abun)

        
        x = self.zero_filled_u_l[0]
        y = self.zero_filled_u_l[1]

        # get the Fourier transform of x and y
        x_dft = rfft(x).real
        y_dft = rfft(y).real

        s_dft_xy = dot(x_dft, y_dft)/(norm(x_dft)*norm(y_dft))

        # using the existing weighted_cosine_correlation function to get S_WC(X,Y)
        s_wc_x_y = self.weighted_cosine_correlation()

        # final step
        s_dft = (n_x * s_wc_x_y + self.n_x_y * s_dft_xy) / (n_x + self.n_x_y)

        return s_dft

    def jaccard_distance(self):
        
        def jaccard_similarity(list1, list2):
            intersection = len(list(set(list1).intersection(list2)))
            union = (len(list1) + len(list2)) - intersection
            return float(intersection) / union
        
       
        correlation = jaccard_similarity(self.zero_filled_u_l[0], self.zero_filled_u_l[1])
        return correlation

    def dwt_correlation(self):
        
        if self.n_x_y == 0: return 0

        # count number of non-zero abundance/peak intensity values
        n_x = sum(a != 0 for a in self.exp_abun)

        #calculate cosine correlation, 
        x = self.zero_filled_u_l[0]
        y = self.zero_filled_u_l[1]

        #Make x and y into an array
        x_a=list(x)
        y_a=list(y)

        # get the wavelet transform of x and y (Daubechies with a filter length of 4. Asymmetric. pywavelets function)
        #Will only use the detail dwt (dwtDd
        (x_dwtA,x_dwtD) = dwt(x_a,'db2')
        (y_dwtA,y_dwtD) = dwt(y_a,'db2')

        s_dwt_xy = dot(x_dwtD, y_dwtD)/(norm(x_dwtD)*norm(y_dwtD))

        # using the existing weighted_cosine_correlation function to get S_WC(X,Y)
        s_wc_x_y = self.weighted_cosine_correlation()

        # final step
        s_dwt = (n_x * s_wc_x_y + self.n_x_y * s_dwt_xy) / (n_x + self.n_x_y)

        return s_dwt