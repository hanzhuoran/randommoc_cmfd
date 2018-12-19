# randommoc_cmfd
A toy model of applying CMFD on random MOC
# How to run code
`python latticerr.py`
# XS data
Data are stored in homo_XS and XS
For now, I have 1 group and 2 group data for homogenized slab. 
For heterogeneos fuel pin and moderator, I have 2 group and 10 group.
# Change paramerter
All the parameters are in util.py
Names of all parameters are very straight forward. 
`DZ` is Dead Zone.
`MAX_D` is the termination distance for a ray.
`xmbd` means `x minus boundary`, which is the left boundary. 
`transreflect` means this boundary is transmissive, but the ray crossing the surface will be bounced back and assign zero angular flux.
`ratio` is to condense multiple cells into an even bigger cell. This functionality is not well tested yet.
