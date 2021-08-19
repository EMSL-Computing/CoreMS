
class KendrickGrouping:

    def mz_odd_even_index_lists(self):
        # get odd and even indexes lists
        
        even_idx = [] 
        odd_idx = []
        
        for i, mspeak in enumerate(self.mspeaks):
            
            if mspeak.nominal_mz_exp % 2 == 0:
                even_idx.append(i)
            else:
                odd_idx.append(i)
        
        return even_idx, odd_idx 

    def calc_error(self, current, test):
        
        return ((test-current)/current)*1e6
  
    
    def populate_kendrick_index_dict_error(self, list_indexes, sort=True):
        
        def error():

            return  abs(current_kmd_reference - next_mspeak.kmd)

        already_found = []
        
        all_results = []
        
        for i in list_indexes:
            
            result_indexes = []
            
            mspeak = self.mspeaks[i]    

            current_kmd_reference = mspeak.kmd
            for j in list_indexes :
                
                if j not in already_found and j != i:
                    
                    next_mspeak = self.mspeaks[j]

                    if  error() <= 0.001:
                        
                        result_indexes.append(j)    
                        already_found.append(j)

                        current_kmd_reference = next_mspeak.kmd
            
            if result_indexes and len(result_indexes) > 3:
                
                already_found.append(i)
                
                result_indexes.insert(0,i)
                
                all_results.append(result_indexes)        
            else:
                
                for w in result_indexes:

                    already_found.remove(w)        

        kendrick_group_index = { i : indexes_list for i, indexes_list in enumerate(all_results) }

       
            #return dictionary with the keys sorted by sum of the abundances
        if sort:
            
            return dict(sorted(kendrick_group_index.items(), key = lambda it: sum([self.mspeaks[i].abundance for i in it[1]]), reverse=False ))
        
        else:
            
            return kendrick_group_index
        
    def populate_kendrick_index_dict_rounding(self, list_indexes, sort=True):

        #self.test(list_indexes)
        #breakpoint()
        kendrick_group_index = {}
        
        for i in list_indexes:
            
            mspeak = self.mspeaks[i]
            
            group = round(mspeak.kmd * 100)
            
            if group not in kendrick_group_index:
                
                kendrick_group_index[group] = [i]

            else: 
                
                last_index = kendrick_group_index[group][-1]
                
                print(abs(mspeak.kmd - self.mspeaks[last_index].kmd ))
                
                if abs(mspeak.kmd - self.mspeaks[last_index].kmd ) < 0.001:

                    kendrick_group_index[group].append(i)
                
               


            #return dictionary with the keys sorted by sum of the abundances
        if sort:
            return dict(sorted(kendrick_group_index.items(), key = lambda it: sum([self.mspeaks[i].abundance for i in it[1]]), reverse=True ))
        
        else:
            return kendrick_group_index

    def sort_abundance_kendrick_dict(self, even_kendrick_group_index, odd_kendrick_group_index):
    
        all_even_indexes = [i for v in even_kendrick_group_index.values() for i in v]
        
        all_odd_indexes = [i for v in odd_kendrick_group_index.values() for i in v]

        sum_even = sum([self.mspeaks[i].abundance for i in all_even_indexes])
        
        sum_odd = sum([self.mspeaks[i].abundance for i in all_odd_indexes])
        
        if sum_even >= sum_odd:
            
            even_kendrick_group_index.update(odd_kendrick_group_index)
            
            return even_kendrick_group_index
        
        else: 
            
            odd_kendrick_group_index.update(even_kendrick_group_index)
            
            return odd_kendrick_group_index  
        
    
    def kendrick_groups_indexes(self, sort=True):
        
        #return dictionary with the kmd as keys and the correspondents peaks indexes
        even_idx, odd_idx = self.mz_odd_even_index_lists()
        
        even_kendrick_group_index = self.populate_kendrick_index_dict_error(even_idx, sort=sort)
        
        odd_kendrick_group_index = self.populate_kendrick_index_dict_error(odd_idx, sort=sort)
        
        return self.sort_abundance_kendrick_dict(even_kendrick_group_index, odd_kendrick_group_index)