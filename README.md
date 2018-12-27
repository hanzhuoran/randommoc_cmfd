# CMFD on Random MOC
This is a course project for MIT reactor physics 2: 22.212
The project is a simplified prototype of applying CMFD on random MOC.
The course is finished, so this repo is currently in libo. 
I will come back on this project and work on it as a research topic after Feburary.
# How to run code
`python latticerr.py`
Set `CMFD` value to 1 to turn it on. (Only work for reflective boundary cases now)
Set `plot` value to 1 to see the trajectories of all rays. (Note: It can be super SLOW!)
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
# Documentation
Methodology and implementation process are documented in the course final report.
Further updates on this project will be posted once I make enough progress.
