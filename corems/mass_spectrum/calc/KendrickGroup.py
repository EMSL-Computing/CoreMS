
class KendrickGrouping:

    def odd_even_index_lists(self):
        # get odd and even indexes lists
        
        even_idx = [] 
        odd_idx = []
        
        for i, mspeak in enumerate(self.mspeaks):
            
            if mspeak.mz_exp % 2 == 0:
                even_idx.append(i)
            else:
                odd_idx.append(i)
        
        return even_idx, odd_idx 

    def populate_kendrick_index_dict(self, list_indexes, sort=True):

        kendrick_group_index = {}
        #set mass spectrum to iterate only over odd numbers
        self.set_indexes(list_indexes)  
        
        for i, mspeak in enumerate(self.mspeaks):
            
            if mspeak.kmd not in kendrick_group_index:
                
                kendrick_group_index[mspeak.kmd] = [i]

            else: 

                kendrick_group_index[mspeak.kmd].append(i)

        self.reset_indexes()
        
            #return dictionary with the keys sorted by sum of the abundances
        if sort:
            
            return sorted(kendrick_group_index.items(), key = lambda it: sum([self.mspeaks[i].abundance for i in it[1]]), reverse=True )
        
        else:
            return kendrick_group_index

    def sort_abundance_kendrick_dict(self, even_kendrick_group_index, odd_kendrick_group_index):
    
        all_even_indexes = [i for v in even_kendrick_group_index.values() for i in v]
        
        all_odd_indexes = [i for v in even_kendrick_group_index.values() for i in v]

        sum_even = sum([self.mspeaks[i].abundance for i in all_even_indexes])
        
        sum_odd = sum([self.mspeaks[i].abundance for i in all_odd_indexes])
        
        if sum_even >= sum_odd:
            
            return even_kendrick_group_index, odd_kendrick_group_index
        
        else: 
            
            return odd_kendrick_group_index, even_kendrick_group_index    
        
    
    def kendrick_groups_indexes(self, sort=True):
        
        #return dictionary with the kmd as keys and the correspondents peaks indexes
        even_idx, odd_idx = self.odd_even_index_lists()

        even_kendrick_group_index = self.populate_kendrick_index_dict(even_idx, sort=sort)
        
        odd_kendrick_group_index = self.populate_kendrick_index_dict(odd_idx, sort=sort)
        
        return self.sort_abundance_kendrick_dict(even_kendrick_group_index, odd_kendrick_group_index)