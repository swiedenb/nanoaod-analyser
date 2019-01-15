def read_in_arrays(datapath):
    tree = uproot.open(datapath)["Events"]
    tree_dict = {}
    for key in tree.keys():
        tree_dict[key] = tree.array(key)
    return tree_dict
    
