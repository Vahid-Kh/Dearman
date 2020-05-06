import scipy.io as spio
import matplotlib.pyplot as plt

mat = spio.loadmat('OD1_prognosis.mat', squeeze_me=True)

"""X"""
OD1CV = mat['OD_1_Count_vect']
"""Y"""
OD1CT = mat['OD_1_Calculation_time']*60*7*24
OD1ST = mat['OD_1_Serv_Time_estimate']*60*7*24

OD1RV = mat['OD_1_Reliability_vect']
OD1FV = mat['OD_1_filt_vect']
label = ['OD_1_Count_vect', 'OD_1_Calculation_time', 'OD_1_Serv_Time_estimate', 'OD_1_Reliability_vect', 'OD_1_filt_vect']
"""___________________________________________________________________________________________________________
____________________________________________________PLOTS_____________________________________________________
___________________________________________________________________________________________________________"""
""" 
First plots the lines and then 'ro' spesifies red dots as the TDN states __ r ::  for red and o :: for circle
    Color               Shape                  shape
    b : blue            "8"	: octagon          "," : pixel 
    g : green           "s"	: square           "o" : circle 
    r : red             "p"	: pentagon         "v" : triangle_down 
    c : cyan            "P"	: plus (filled)    "x" : x
    m : magenta         "*"	: star             "X" : x (filled)
    y : yellow          "h"	: hexagon1         "D" : diamond
    k : black           "H"	: hexagon2         "d" : thin_diamond  
    w : white           "+"	: plus                       

"""

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

fig, host = plt.subplots()
fig.subplots_adjust(right=0.75)

par1 = host.twinx()
par2 = host.twinx()
par3 = host.twinx()

# Offset the right spine of par2.  The ticks and label have already been
# placed on the right by twinx above.
par1.spines["right"].set_position(("axes", 0))
par2.spines["right"].set_position(("axes", 1))

# Having been created by twinx, par2 has its frame off, so the line of its
# detached spine is invisible.  First, activate the frame but make the patch
# and spines invisible.
make_patch_spines_invisible(par2)
par3.spines["right"].set_position(("axes", 1.1))
make_patch_spines_invisible(par3)



# Second, show the right spine.
par1.spines["right"].set_visible(True)
par2.spines["right"].set_visible(True)
par3.spines["right"].set_visible(True)


p1, = host.plot(OD1CV, OD1CT, "b-", label=label[1])
p2, = par1.plot(OD1CV, OD1ST, "r-", label=label[2])
plt.grid(True)
p3, = par2.plot(OD1CV, OD1RV, "g-", label=label[3])
p4, = par3.plot(OD1CV, OD1FV, "c-", label=label[4])

host.set_xlim(min(OD1CV), max(OD1CV))
host.set_ylim(min(min(OD1CT), min(OD1ST)) * 0.9, max(max(OD1CT), max(OD1ST)) * 1.1)
par1.set_ylim(min(min(OD1CT), min(OD1ST)) * 0.9, max(max(OD1CT), max(OD1ST)) * 1.1)
par2.set_ylim(min(OD1RV) * 0.9, max(OD1RV) * 1.1)
par3.set_ylim(min(OD1FV) * 0.85, max(OD1FV) * 1.15)

host.set_xlabel(label[0])
host.set_ylabel(label[1])
par1.set_ylabel(label[2])
par2.set_ylabel(label[3])
par3.set_ylabel(label[4])

host.yaxis.label.set_color(p1.get_color())
par1.yaxis.label.set_color(p2.get_color())
par2.yaxis.label.set_color(p3.get_color())
par3.yaxis.label.set_color(p4.get_color())

tkw = dict(size=4, width=1.5)
host.tick_params(axis='y', colors=p1.get_color(), **tkw)
par1.tick_params(axis='y', colors=p2.get_color(), **tkw)
par2.tick_params(axis='y', colors=p3.get_color(), **tkw)
par3.tick_params(axis='y', colors=p4.get_color(), **tkw)
host.tick_params(axis='x', **tkw)

lines = [p1, p2, p3, p4]

host.legend(lines, [l.get_label() for l in lines])
plt.title('System reliability as a function of time')

plt.show()