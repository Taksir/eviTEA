from stages import *
rp=pickle.load(open(WORK/"g1d.pkl","rb"))
fill_dollars(sorted(rp.selected_sequences), budget=155)
